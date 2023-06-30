#include <elaspect/material_model/handler.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  template <int dim>
  void MaterialHandler<dim>::initialize()
  {
    basic_property_model->initialize();

    if (viscosity_model.get())
      viscosity_model->initialize();

    if (plasticity_model.get())
      plasticity_model->initialize();
  }


  template <int dim>
  void MaterialHandler<dim>::update()
  {
    basic_property_model->update();

    if (viscosity_model.get())
      viscosity_model->update();

    if (plasticity_model.get())
      plasticity_model->update();
  }


  template <int dim>
  void
  MaterialHandler<dim>::
  evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
           MaterialModel::MaterialModelOutputs<dim> &out) const
  {
    const unsigned int n_points = in.n_evaluation_points();
    const unsigned int n_fields = this->n_compositional_fields() + 1;

    MaterialModel::BasicPropertyOutputs<dim> basic_property_outputs(n_fields);
    std::vector<double> volume_fractions(n_fields);
    std::vector<double> K(n_fields), Gve(n_fields), eta(n_fields);
    SymmetricTensor<4,dim> D_avg;
    SymmetricTensor<2,dim> sigma_old, s_old, depsilon, sigma_tr;

    for (unsigned int i = 0; i < n_points; ++i)
    {
      volume_fractions = MaterialUtilities::compute_composition_fractions(in.composition[i]);

      // calculate basic properties if requested
      if (in.requested_properties & MaterialModel::MaterialProperties::basic_properties)
      {
        basic_property_model->compute_basic_properties(in, i, basic_property_outputs);

        if (in.requested_properties & MaterialModel::MaterialProperties::specific_heat)
          out.specific_heat[i] = MaterialUtilities::average_value(volume_fractions, 
                                                                  basic_property_outputs.specific_heat_capacities, 
                                                                  MaterialUtilities::arithmetic);
        if (in.requested_properties & MaterialModel::MaterialProperties::thermal_expansivity)
          out.thermal_expansivities[i] = MaterialUtilities::average_value(volume_fractions, 
                                                                          basic_property_outputs.thermal_expansivities, 
                                                                          MaterialUtilities::arithmetic);
        if (in.requested_properties & MaterialModel::MaterialProperties::thermal_conductivity)
          out.thermal_conductivities[i] = MaterialUtilities::average_value(volume_fractions, 
                                                                           basic_property_outputs.thermal_conductivities, 
                                                                           MaterialUtilities::arithmetic);
        if (in.requested_properties & MaterialModel::MaterialProperties::density)
          out.densities[i] = MaterialUtilities::average_value(volume_fractions, 
                                                              basic_property_outputs.densities, 
                                                              MaterialUtilities::arithmetic);
        if (in.requested_properties & MaterialModel::MaterialProperties::reaction_terms)
          out.reaction_terms[i] = basic_property_outputs.reaction_terms;
      }

      if (!(in.requested_properties & MaterialModel::MaterialProperties::tangent_modulus))
      {
        if (in.requested_properties & MaterialModel::MaterialProperties::viscosity)
        {
          for (unsigned int field = 0; field < n_fields; ++field)
            eta[field] = viscosity_model->compute_viscosity(in, i, field);
          out.viscosities[i] = MaterialUtilities::average_value (volume_fractions, eta, rheological_parameters_averaging);
        }

        continue;
      }

      for (unsigned int field = 0; field < n_fields; ++field)
      {
        // calculate the stiffness modulus for each field
        K[field]   = basic_property_outputs.bulk_moduli[field];
        Gve[field] = basic_property_outputs.shear_moduli[field];
        if (this->constitutive_relation() & ConstitutiveRelation::viscosity)
        {
          eta[field] = viscosity_model->compute_viscosity(in, i, field);
          Gve[field] = (Gve[field] * eta[field]) / (Gve[field] * this->get_timestep() + eta[field]);
        }
      }

      // do averaging for viscosity and elastic moduli
      if (in.requested_properties & MaterialModel::MaterialProperties::viscosity)
        out.viscosities[i] = MaterialUtilities::average_value(volume_fractions, eta, 
                                                              rheological_parameters_averaging);

      const double K_avg   = MaterialUtilities::average_value(volume_fractions, K, 
                                                            rheological_parameters_averaging);
      const double Gve_avg = MaterialUtilities::average_value(volume_fractions, Gve, 
                                                              rheological_parameters_averaging);
      MaterialUtilities::compute_isotropic_stiffness_tensor(K_avg, Gve_avg, D_avg);

      if (this->constitutive_relation() & ConstitutiveRelation::plasticity)
      {
        // Return-mapping cannot be executed on faces or subfaces.
        AssertDimension(n_points, this->get_qpd_handler().n_quadrature_points());

        // calculate the consistent tangent modulus
        sigma_old = in.old_stress[i];
        depsilon  = symmetrize(in.incremental_displacement_gradient[i]);

        if (this->constitutive_relation() & ConstitutiveRelation::viscosity)
        {
          const double G_avg = MaterialUtilities::average_value(volume_fractions,
                                                                basic_property_outputs.shear_moduli,
                                                                rheological_parameters_averaging);

          const double p_old = trace(sigma_old) / dim;
          s_old = deviator(sigma_old) * (Gve_avg / G_avg);
          sigma_old = s_old + p_old * unit_symmetric_tensor<dim>();
        }

        if (this->constitutive_relation() & ConstitutiveRelation::thermal_expansion)
        {
          const double alpha = MaterialUtilities::average_value(volume_fractions,
                                                                basic_property_outputs.thermal_expansivities,
                                                                MaterialUtilities::arithmetic);
          const double T_old = in.qpd_cell->get_scalar(i, this->introspection().qpd_indicators.old_temperature);
          const double alpha_dT = alpha * (in.temperature[i] - T_old);
          for (unsigned int d = 0; d < dim; ++d)
            depsilon[d][d] -= alpha_dT;
        }

        sigma_tr = sigma_old + D_avg * depsilon;

        plasticity_model->compute_tangent_modulus(in, i, volume_fractions, 
                                                  sigma_tr, D_avg, out.tangent_moduli[i]);
      }
      else
        out.tangent_moduli[i] = D_avg;
    }
  }


  template <int dim>
  void
  MaterialHandler<dim>::
  apply_return_mapping(const MaterialModel::MaterialModelInputs<dim> &in,
                       std::vector<SymmetricTensor<2,dim>> &          stresses,
                       std::vector<double> &                          plastic_strains,
                       std::vector<double> &                          flags) const
  {
    Assert(this->constitutive_relation() & ConstitutiveRelation::plasticity,
           ExcMessage("MaterialModel::Interface<dim>::execute_return_mapping() "
                      "should not be called when plasticity is not included in the "
                      "constitutive relations."));

    const unsigned int n_fields = this->n_compositional_fields() + 1;
    const unsigned int n_points = in.n_evaluation_points();

    // Return-mapping cannot be executed on faces or subfaces.
    AssertDimension(n_points, this->get_qpd_handler().n_quadrature_points());

    MaterialModel::BasicPropertyOutputs<dim> basic_property_outputs(n_fields);
    std::vector<double> volume_fractions(n_fields);
    std::vector<double> K(n_fields), Gve(n_fields), eta(n_fields);
    SymmetricTensor<4,dim> D_avg;
    SymmetricTensor<2,dim> sigma_old, s_old, depsilon, sigma_tr;
    MaterialModel::Plasticity::StateUpdates<dim> state_updates;

    for (unsigned int i = 0; i < n_points; ++i)
    {
      volume_fractions = MaterialUtilities::compute_composition_fractions(in.composition[i]);

      // calculate the (visco)elastic modulus
      basic_property_model->compute_basic_properties(in, i, basic_property_outputs);
      for (unsigned int field = 0; field < n_fields; ++field)
      {
        K[field]   = basic_property_outputs.bulk_moduli[field];
        Gve[field] = basic_property_outputs.shear_moduli[field];
        if (this->constitutive_relation() & ConstitutiveRelation::viscosity)
        {
          eta[field] = viscosity_model->compute_viscosity(in, i, field);
          Gve[field] = (Gve[field] * eta[field]) / (Gve[field] * this->get_timestep() + eta[field]);
        }
      }

      // do averaging for elastic moduli
      const double K_avg   = MaterialUtilities::average_value(volume_fractions, K, 
                                                              rheological_parameters_averaging);
      const double Gve_avg = MaterialUtilities::average_value(volume_fractions, Gve, 
                                                              rheological_parameters_averaging);
      MaterialUtilities::compute_isotropic_stiffness_tensor(K_avg, Gve_avg, D_avg);

      // calculate trial stress
      sigma_old = in.old_stress[i];
      depsilon  = symmetrize(in.incremental_displacement_gradient[i]);

      if (this->constitutive_relation() & ConstitutiveRelation::viscosity)
      {
        const double G_avg = MaterialUtilities::average_value(volume_fractions,
                                                              basic_property_outputs.shear_moduli,
                                                              rheological_parameters_averaging);

        const double p_old = trace(sigma_old) / dim;
        s_old = deviator(sigma_old) * (Gve_avg / G_avg);
        sigma_old = s_old + p_old * unit_symmetric_tensor<dim>();
      }

      if (this->constitutive_relation() & ConstitutiveRelation::thermal_expansion)
      {
        const double alpha = MaterialUtilities::average_value(volume_fractions, 
                                                              basic_property_outputs.thermal_expansivities, 
                                                              MaterialUtilities::arithmetic);
        const double T_old = in.qpd_cell->get_scalar(i, this->introspection().qpd_indicators.old_temperature);
        const double alpha_dT = alpha * (in.temperature[i] - T_old);
        for (unsigned int d = 0; d < dim; ++d)
          depsilon[d][d] -= alpha_dT;
      }

      sigma_tr = sigma_old + D_avg * depsilon;

      plasticity_model->apply_return_mapping(in, i, volume_fractions,
                                             sigma_tr, D_avg, state_updates);

      stresses[i]        = state_updates.stress;
      plastic_strains[i] = state_updates.plastic_strain;
      flags[i]           = state_updates.flag;
    }
  }


  template <int dim>
  MaterialModel::FieldDependences::Dependence
  MaterialHandler<dim>::
  get_field_dependences_for_evaluation(const MaterialModel::MaterialProperties::Property properties) const
  {
    MaterialModel::FieldDependences::Dependence 
    dependences = MaterialModel::FieldDependences::composition;

    if (properties & MaterialModel::MaterialProperties::basic_properties)
      dependences |= basic_property_model->get_field_dependences(properties);

    if (!(properties & MaterialModel::MaterialProperties::tangent_modulus))
    {
      if (properties & MaterialModel::MaterialProperties::viscosity)
      {
        AssertThrow (this->constitutive_relation() & ConstitutiveRelation::viscosity,
                     ExcMessage("Material property 'viscosity' is requested by the material model inputs, "
                                "but viscosity is not included in the constitutive relations."));

        dependences |= viscosity_model->get_field_dependences();
      }
    }
    else
    {
      if (this->constitutive_relation() & ConstitutiveRelation::thermal_expansion)
        dependences |= MaterialModel::FieldDependences::temperature;

      if (this->constitutive_relation() & ConstitutiveRelation::viscosity)
        dependences |= viscosity_model->get_field_dependences();

      if (this->constitutive_relation() & ConstitutiveRelation::plasticity)
      {
        dependences |= plasticity_model->get_field_dependences();
        dependences |= MaterialModel::FieldDependences::displacement_gradient |
                       MaterialModel::FieldDependences::old_stress;
      }
    }

    return dependences;
  }


  template <int dim>
  MaterialModel::FieldDependences::Dependence
  MaterialHandler<dim>::get_field_dependences_for_return_mapping() const
  {
    AssertThrow (this->constitutive_relation() & ConstitutiveRelation::plasticity,
                 ExcInternalError());

    MaterialModel::FieldDependences::Dependence 
    dependences = MaterialModel::FieldDependences::displacement_gradient |
                  MaterialModel::FieldDependences::composition |
                  MaterialModel::FieldDependences::old_stress;

    const MaterialModel::MaterialProperties::Property 
    properties = MaterialModel::MaterialProperties::tangent_modulus;
    dependences |= basic_property_model->get_field_dependences(properties);

    if (this->constitutive_relation() & ConstitutiveRelation::viscosity)
      dependences |= viscosity_model->get_field_dependences();

    if (this->constitutive_relation() & ConstitutiveRelation::thermal_expansion)
      dependences |= MaterialModel::FieldDependences::temperature;

    return dependences;
  }


  template <int dim>
  MaterialModel::MaterialProperties::Property
  MaterialHandler<dim>::get_material_properties_for_return_mapping () const
  {
    AssertThrow (this->constitutive_relation() & ConstitutiveRelation::plasticity,
                 ExcInternalError());

    MaterialModel::MaterialProperties::Property 
    properties = MaterialModel::MaterialProperties::tangent_modulus;

    if (this->constitutive_relation() & ConstitutiveRelation::viscosity)
      properties |= MaterialModel::MaterialProperties::viscosity;

    if (this->constitutive_relation() & ConstitutiveRelation::thermal_expansion)
      properties |= MaterialModel::MaterialProperties::thermal_expansivity;

    return properties;
  }


  template <int dim>
  void
  MaterialHandler<dim>::declare_parameters (ParameterHandler &prm)
  {
    prm.enter_subsection ("Material model");
    {
      prm.declare_entry ("Rheological parameters averaging scheme", "maximum composition",
                         Patterns::Selection("arithmetic|harmonic|geometric|maximum composition"),
                         "");
    }
    prm.leave_subsection ();

    MaterialModel::BasicProperties::declare_parameters<dim> (prm);
    MaterialModel::Viscosity::declare_parameters<dim> (prm);
    MaterialModel::Plasticity::declare_parameters<dim> (prm);
  }


  template <int dim>
  void
  MaterialHandler<dim>::parse_parameters (ParameterHandler &prm)
  {
    prm.enter_subsection ("Material model");
    {
      rheological_parameters_averaging =
        MaterialUtilities::parse_compositional_averaging_operation (
          "Rheological parameters averaging scheme", prm);
    }
    prm.leave_subsection ();

    basic_property_model.reset (MaterialModel::BasicProperties::create_basic_property_model<dim>(prm));
    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(basic_property_model.get()))
      sim->initialize_simulator (this->get_simulator());
    basic_property_model->parse_parameters (prm);
    basic_property_model->initialize ();

    if (this->constitutive_relation() & ConstitutiveRelation::viscosity)
    {
      viscosity_model.reset (MaterialModel::Viscosity::create_viscosity_model<dim>(prm));
      if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(viscosity_model.get()))
        sim->initialize_simulator (this->get_simulator());
      viscosity_model->parse_parameters (prm);
      viscosity_model->initialize ();
    }

    if (this->constitutive_relation() & ConstitutiveRelation::plasticity)
    {
      plasticity_model.reset (MaterialModel::Plasticity::create_plasticity_model<dim>(prm));
      if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(plasticity_model.get()))
        sim->initialize_simulator (this->get_simulator());
      plasticity_model->parse_parameters (prm);
      plasticity_model->initialize ();
    }
  }
}


// explicit instantiations
namespace elaspect
{
#define INSTANTIATE(dim) \
  template class MaterialHandler<dim>; 

  ELASPECT_INSTANTIATE(INSTANTIATE)
}
