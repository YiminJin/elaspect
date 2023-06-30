#include <elaspect/material_model/io_interface.h>
#include <elaspect/introspection.h>

namespace elaspect
{
  namespace MaterialModel
  {
    // --------------------- MaterialModelInputs ----------------------

    template <int dim>
    MaterialModelInputs<dim>::MaterialModelInputs(const unsigned int                  n_points,
                                                  const unsigned int                  n_comp,
                                                  const FieldDependences::Dependence  dependences,
                                                  const MaterialProperties::Property  properties)
      :
      position(dependences & FieldDependences::position ? n_points : 0, 
               Point<dim>(numbers::signaling_nan<Tensor<1,dim>>())),
      incremental_displacement(dependences & FieldDependences::displacement ? n_points : 0,
                               numbers::signaling_nan<Tensor<1,dim>>()),
      incremental_displacement_gradient(dependences & FieldDependences::displacement_gradient ? n_points : 0, 
                                        numbers::signaling_nan<Tensor<2,dim>>()),
      fluid_pressure(dependences & FieldDependences::fluid_pressure ? n_points : 0,
                     numbers::signaling_nan<double>()),
      temperature(dependences & FieldDependences::temperature ? n_points : 0, 
                  numbers::signaling_nan<double>()),
      composition(dependences & FieldDependences::composition ? n_points : 0,
                  std::vector<double>(n_comp, numbers::signaling_nan<double>())),
      old_stress(dependences & FieldDependences::old_stress ? n_points : 0,
                 numbers::signaling_nan<SymmetricTensor<2,dim>>()),
      qpd_cell(),
      field_dependences(dependences),
      requested_properties(properties)
    {
      AssertThrow(field_dependences != FieldDependences::none,
                  ExcMessage("The MaterialModelInputs object is empty!"));
    }


    template <int dim>
    MaterialModelInputs<dim>::
    MaterialModelInputs(const DataPostprocessorInputs::Vector<dim> &input_data,
                        const Introspection<dim>                   &introspection,
                        const FieldDependences::Dependence          dependences,
                        const MaterialProperties::Property          properties)
      :
      position(input_data.evaluation_points),
      incremental_displacement(dependences & FieldDependences::displacement ? 
                               input_data.solution_values.size() : 0,
                               numbers::signaling_nan<Tensor<1,dim>>()),
      incremental_displacement_gradient(dependences & FieldDependences::displacement_gradient ?
                                        input_data.solution_values.size() : 0,
                                        numbers::signaling_nan<Tensor<2,dim>>()),
      fluid_pressure(dependences & FieldDependences::fluid_pressure ? 
                     input_data.solution_values.size() : 0,
                     numbers::signaling_nan<double>()),
      temperature(dependences & FieldDependences::temperature ?
                  input_data.solution_values.size() : 0, 
                  numbers::signaling_nan<double>()),
      composition(dependences & FieldDependences::composition ?
                  input_data.solution_values.size() : 0,
                  std::vector<double>(introspection.n_compositional_fields,
                                      numbers::signaling_nan<double>())),
      old_stress(dependences & FieldDependences::old_stress ?
                 input_data.solution_values.size(): 0,
                 numbers::signaling_nan<SymmetricTensor<2,dim>>()),
      qpd_cell(),
      field_dependences(dependences),
      requested_properties(properties)
    {
      const unsigned int n_points = input_data.solution_values.size();

      if (field_dependences & FieldDependences::displacement)
      {
        for (unsigned int q = 0; q < n_points; ++q)
          for (unsigned int d = 0; d < dim; ++d)
            incremental_displacement[q][d] = input_data.solution_values[q][introspection.component_indices.displacement[d]];
      }

      if (field_dependences & FieldDependences::displacement_gradient)
      {
        for (unsigned int q = 0; q < n_points; ++q)
        {
          Tensor<2,dim> grad_u;
          for (unsigned int d = 0; d < dim; ++d)
          {
            grad_u[d] = input_data.solution_gradients[q][introspection.component_indices.displacement[d]];
            incremental_displacement_gradient[q] = symmetrize(grad_u);
          }
        }
      }

      if (field_dependences & FieldDependences::temperature)
      {
        for (unsigned int q = 0; q < n_points; ++q)
          temperature[q] = input_data.solution_values[q][introspection.component_indices.temperature];
      }

      if (field_dependences & FieldDependences::fluid_pressure)
      {
        const unsigned int fluid_pressure_index = introspection.variable("fluid pressure").first_component_index;
        for (unsigned int q = 0; q < n_points; ++q)
          fluid_pressure[q] = input_data.solution_values[q][fluid_pressure_index];
      }

      if (field_dependences & FieldDependences::composition)
      {
        for (unsigned int q = 0; q < n_points; ++q)
          for (unsigned int i = 0; i < introspection.n_compositional_fields; ++i)
            composition[q][i] = input_data.solution_values[q][introspection.component_indices.compositional_fields[i]];
      }

      if (field_dependences & FieldDependences::old_stress)
      {
        // When doing postprocessing, we use the current stress value to compute
        // material properties that depend on stress, such as density and viscosity.
        for (unsigned int q = 0; q < n_points; ++q)
          for (unsigned int c = 0; c < SymmetricTensor<2,dim>::n_independent_components; ++c)
            old_stress[q].access_raw_entry(c) = input_data.solution_values[q][introspection.component_indices.stress[c]];
      }
    }


    template <int dim>
    void
    MaterialModelInputs<dim>::reinit (const FEValuesBase<dim>                   &fe_values,
                                      const QPDHandler<dim>                     &qpd_handler,
                                      const Introspection<dim>                  &introspection,
                                      const TrilinosWrappers::MPI::BlockVector  &solution)
    {
      if (field_dependences & FieldDependences::position)
        position = fe_values.get_quadrature_points();

      if (field_dependences & FieldDependences::displacement)
        fe_values[introspection.extractors.displacement].get_function_values(solution, incremental_displacement);

      if (field_dependences & FieldDependences::displacement_gradient)
        fe_values[introspection.extractors.displacement].get_function_gradients(solution, incremental_displacement_gradient);

      if (field_dependences & FieldDependences::temperature)
        fe_values[introspection.extractors.temperature].get_function_values(solution, temperature);

      if (field_dependences & FieldDependences::fluid_pressure)
      {
        const FEValuesExtractors::Scalar fp_extractor(introspection.variable("fluid pressure").first_component_index);
        fe_values[fp_extractor].get_function_values (solution, fluid_pressure);
      }

      if (dynamic_cast<const FEValues<dim>*>(&fe_values) != nullptr)
      {
        AssertDimension(fe_values.n_quadrature_points, qpd_handler.n_quadrature_points());

        // get quadrature point data directly from QPDHandler
        qpd_cell = typename QPDHandler<dim>::active_cell_iterator(*fe_values.get_cell(),
                                                                  &qpd_handler);

        if (field_dependences & FieldDependences::composition)
        {
          std::vector<double> component_values(fe_values.n_quadrature_points);
          for (unsigned int i = 0; i < introspection.n_compositional_fields; ++i)
          {
            qpd_cell->get(introspection.qpd_indicators.composition + i, component_values);
            for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q)
              composition[q][i] = component_values[q];
          }
        }

        if (field_dependences & FieldDependences::old_stress)
        {
          for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q)
            old_stress[q] = qpd_cell->get_symmetric_tensor(q, introspection.qpd_indicators.old_stress);
        }
      }
      else
      {
        // get quadrature point data from the finite element space
        if (field_dependences & FieldDependences::composition)
        {
          std::vector<double> component_values(fe_values.n_quadrature_points);
          for (unsigned int i = 0; i < introspection.n_compositional_fields; ++i)
          {
            fe_values[introspection.extractors.compositional_fields[i]].get_function_values(solution, component_values);
            for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q)
              composition[q][i] = component_values[q];
          }
        }

        if (field_dependences & FieldDependences::old_stress)
          fe_values[introspection.extractors.stress].get_function_values(solution, old_stress);
      }
    }


    template <int dim>
    unsigned int
    MaterialModelInputs<dim>::n_evaluation_points() const
    {
      if (field_dependences & FieldDependences::position)
        return position.size();

      if (field_dependences & FieldDependences::displacement)
        return incremental_displacement.size();

      if (field_dependences & FieldDependences::displacement_gradient)
        return incremental_displacement_gradient.size();

      if (field_dependences & FieldDependences::fluid_pressure)
        return fluid_pressure.size();

      if (field_dependences & FieldDependences::temperature)
        return temperature.size();

      AssertThrow (false, ExcInternalError());
    }


    // --------------------- MaterialModelOutputs ----------------------

    template <int dim>
    MaterialModelOutputs<dim>::MaterialModelOutputs (const unsigned int n_points,
                                                     const unsigned int n_comp)
      :
      tangent_moduli(n_points, numbers::signaling_nan<SymmetricTensor<4,dim> >()),
      viscosities(n_points, numbers::signaling_nan<double>()),
      densities(n_points, numbers::signaling_nan<double>()),
      specific_heat(n_points, numbers::signaling_nan<double>()),
      thermal_expansivities(n_points, numbers::signaling_nan<double>()),
      thermal_conductivities(n_points, numbers::signaling_nan<double>()),
      entropy_derivative_temperature(n_points, numbers::signaling_nan<double>()),
      entropy_derivative_pressure(n_points, numbers::signaling_nan<double>()),
      reaction_terms(n_points, std::vector<double>(n_comp, numbers::signaling_nan<double>()))
    {}


    template <int dim>
    MaterialModelOutputs<dim>::
    MaterialModelOutputs (const MaterialModelOutputs<dim> &source)
      :
      tangent_moduli (source.tangent_moduli),
      viscosities (source.viscosities),
      densities (source.densities),
      specific_heat (source.specific_heat),
      thermal_expansivities (source.thermal_expansivities),
      thermal_conductivities (source.thermal_conductivities),
      entropy_derivative_temperature (source.entropy_derivative_temperature),
      entropy_derivative_pressure (source.entropy_derivative_pressure),
      reaction_terms (source.reaction_terms)
    {
      Assert (source.additional_outputs.size() == 0,
              ExcMessage ("You can not copy MaterialModelOutputs objects that have "
                          "additional output objects attached."));
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
    template struct MaterialModelInputs<dim>; \
    template struct MaterialModelOutputs<dim>;

    ELASPECT_INSTANTIATE(INSTANTIATE)
  }
}
