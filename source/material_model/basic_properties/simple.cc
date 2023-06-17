#include <elaspect/material_model/basic_properties/simple.h>
#include <elaspect/utilities.h>

namespace elaspect
{
  namespace MaterialModel
  {
    namespace BasicProperties
    {
      template <int dim>
      void
      Simple<dim>::
      compute_basic_properties (const MaterialModel::MaterialModelInputs<dim> &in,
                                const unsigned int i,
                                MaterialModel::BasicPropertyOutputs<dim> &out) const
      {
        if (in.requested_properties & MaterialProperties::specific_heat)
          out.specific_heat_capacities = specific_heat_capacities;

        if (in.requested_properties & MaterialProperties::thermal_expansivity)
          out.thermal_expansivities = thermal_expansivities;

        if (in.requested_properties & MaterialProperties::thermal_conductivity)
          out.thermal_conductivities = thermal_conductivities;

        if (in.requested_properties & MaterialProperties::tangent_modulus)
        {
          out.bulk_moduli  = bulk_moduli;
          out.shear_moduli = shear_moduli;
        }

        if (in.requested_properties & MaterialProperties::density)
        {
          // rho = rho_ref * (1 - alpha * dim * (T - T_ref) + p / K)
          const double pressure = -1. / dim * trace(in.old_stress[i]);
          const double dT = in.temperature[i] - reference_temperature;
          for (unsigned int field = 0; field < reference_densities.size(); ++field)
            out.densities[field] = reference_densities[field] * 
                                   (1. - thermal_expansivities[field] * dim * dT + pressure / bulk_moduli[field]);
        }

        if (in.requested_properties & MaterialProperties::reaction_terms)
        {
          // no reaction
          for (unsigned int comp = 0; comp < out.reaction_terms.size(); ++comp)
            out.reaction_terms[comp] = 0;
        }
      }


      template <int dim>
      FieldDependences::Dependence
      Simple<dim>::
      get_field_dependences (const MaterialProperties::Property requested_properties) const
      {
        FieldDependences::Dependence dependences = FieldDependences::none;
        if (requested_properties & MaterialProperties::density)
          dependences |= FieldDependences::temperature |
                         FieldDependences::old_stress;;
        return dependences;
      }


      template <int dim>
      void
      Simple<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection ("Material model");
        {
          prm.enter_subsection ("Basic properties");
          {
            prm.enter_subsection ("Simple model");
            {
              prm.declare_entry ("Reference temperature", "273.15",
                                 Patterns::Double(0),
                                 "");

              prm.declare_entry ("Reference densities", "2700",
                                 Patterns::List(Patterns::Double(0)),
                                 "");

              prm.declare_entry ("Heat capacities", "1.25e3",
                                 Patterns::List(Patterns::Double(0)),
                                 "");

              prm.declare_entry ("Thermal expansivities", "1e-5",
                                 Patterns::List(Patterns::Double(0)),
                                 "");

              prm.declare_entry ("Thermal conductivities", "3.0",
                                 Patterns::List(Patterns::Double(0)),
                                 "");

              prm.declare_entry ("Young's moduli", "1.e10",
                                 Patterns::List(Patterns::Double(0)),
                                 "");

              prm.declare_entry ("Poisson's ratio", "0.3",
                                 Patterns::List(Patterns::Double(0)),
                                 "");
            }
            prm.leave_subsection ();
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();
      }


      template <int dim>
      void Simple<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection ("Material model");
        {
          prm.enter_subsection ("Basic properties");
          {
            prm.enter_subsection ("Simple model");
            {
              reference_temperature = prm.get_double("Reference temperature");

              const unsigned int n_fields = this->n_compositional_fields() + 1;

              specific_heat_capacities = Utilities::possibly_extend_from_1_to_N (
                Utilities::string_to_double(Utilities::split_string_list(prm.get("Heat capacities"))),
                n_fields,
                "Heat capacities");

              reference_densities = Utilities::possibly_extend_from_1_to_N (
                Utilities::string_to_double(Utilities::split_string_list(prm.get("Reference densities"))),
                n_fields,
                "Reference densities");

              thermal_expansivities = Utilities::possibly_extend_from_1_to_N (
                Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal expansivities"))),
                n_fields,
                "Thermal expansivities");

              thermal_conductivities = Utilities::possibly_extend_from_1_to_N (
                Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal conductivities"))),
                n_fields,
                "Thermal conductivities");

              std::vector<double> E, nu;
              E = Utilities::possibly_extend_from_1_to_N (
                Utilities::string_to_double(Utilities::split_string_list(prm.get("Young's moduli"))),
                n_fields,
                "Young's moduli");
              nu = Utilities::possibly_extend_from_1_to_N (
                Utilities::string_to_double(Utilities::split_string_list(prm.get("Poisson's ratio"))),
                n_fields,
                "Poisson's moduli");

              bulk_moduli.resize (n_fields);
              shear_moduli.resize (n_fields);
              for (unsigned int i = 0; i < n_fields; ++i)
              {
                bulk_moduli[i] = (dim == 2 ?
                    E[i] / (2 * (1 + nu[i]) * (1 - 2 * nu[i])) :
                    E[i] / (3 * (1 - 2 * nu[i])));

                shear_moduli[i] = E[i] / (2 * (1 + nu[i]));
              }
            }
            prm.leave_subsection ();
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();
      }
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace MaterialModel
  {
    namespace BasicProperties
    {
      ELASPECT_REGISTER_BASIC_PROPERTY_MODEL(Simple,
                                         "simple",
                                         "")
    }
  }
}
