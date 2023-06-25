#include <elaspect/material_model/viscosity/diffusion_dislocation.h>
#include <elaspect/utilities.h>

namespace elaspect
{
  namespace MaterialModel
  {
    namespace Viscosity
    {
      template <int dim>
      double
      DiffusionDislocation<dim>::
      compute_viscosity(const MaterialModel::MaterialModelInputs<dim> &in,
                        const unsigned int i,
                        const unsigned int field) const
      {
        // Power law creep equation
        // edot_ii = A_i * s_ii^n + d^{-m} \exp\left(-\frac{E + pV}{RT}\right)
        // where ii indicates the square root of the second invariant.
        // Thus, the viscosity can be expressed as
        // \eta = s_ii / (2 edot_ii) 
        //      = \frac{1}{2} A_i^{-1} s_ii^{1-n} d^m \exp\left(\frac{E + pV}{RT}\right).
        // Here we use the stress from the previous time step to calculate s_ii and p.
        const double p_old = -trace(in.old_stress[i]) / dim;
        const double s_old = std::sqrt(std::abs(second_invariant(in.old_stress[i])));

        const double RT = constants::gas_constant * in.temperature[i];

        const double big_number = std::sqrt(std::numeric_limits<double>::max());
        double eta_diff = big_number, eta_disl = big_number;
        if (A_diff[field] > 0)
        {
          eta_diff = 0.5 / A_diff[field] * std::pow(d, m_diff[field]) *
                     std::exp((E_diff[field] + p_old * V_diff[field]) / RT);

          AssertThrow(eta_diff > 0,
                      ExcMessage("Negative diffusion viscosity detected. This is unphysical and "
                                 "should not happen. Check for negative parameters. Temperature "
                                 "and pressure are "
                                 + Utilities::to_string(in.temperature[i]) + "K and " 
                                 + Utilities::to_string(p_old) + "Pa."));

          // limit the viscosity to prevent floating-point overflow errors
          eta_diff = std::min(eta_diff, big_number);
        }

        if (A_disl[field] > 0)
        {
          eta_disl = 0.5 / A_disl[field] * std::pow(s_old, 1 - n_disl[field]) *
                     std::exp((E_disl[field] + p_old * V_disl[field]) / RT);

          AssertThrow(eta_disl > 0,
                      ExcMessage("Negative dislocation viscosity detected. This is unphysical and "
                                 "should not happen. Check for negative parameters. Temperature, "
                                 "deviatoric stress and pressure are "
                                 + Utilities::to_string(in.temperature[i]) + "K, " 
                                 + Utilities::to_string(s_old) + "Pa and "
                                 + Utilities::to_string(p_old) + "Pa."));

          // limit the viscosity to prevent floating-point overflow errors
          eta_disl = std::min(eta_disl, big_number);
        }

        return eta_min + 1. / (1. / eta_diff + 1. / eta_disl + 1. / eta_max);
      }


      template <int dim>
      FieldDependences::Dependence
      DiffusionDislocation<dim>::get_field_dependences() const
      {
        return FieldDependences::temperature |
               FieldDependences::old_stress;
      }


      template <int dim>
      void
      DiffusionDislocation<dim>::declare_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Material model");
        {
          prm.enter_subsection("Viscosity");
          {
            prm.enter_subsection("Diffusion dislocation");
            {
              // Parameters for diffusion creep
              prm.declare_entry ("Prefactors for diffusion creep", "1.5e-15",
                                 Patterns::Anything(),
                                 "List of viscosity prefactors, $A$, for background material and compositional fields, "
                                 "for a total of N+1 values, where N is the number of compositional fields. "
                                 "If only one value is given, then all use the same value. "
                                 "Units: \\si{\\per\\pascal\\meter}$^{m_{\\text{diffusion}}}$\\si{\\per\\second}.");
              prm.declare_entry ("Grain size exponents for diffusion creep", "3.",
                                 Patterns::Anything(),
                                 "List of grain size exponents, $m_{\\text{diffusion}}$, for background material and compositional fields, "
                                 "for a total of N+1 values, where N is the number of compositional fields. "
                                 "If only one value is given, then all use the same value. Units: None.");
              prm.declare_entry ("Activation energies for diffusion creep", "375e3",
                                 Patterns::Anything(),
                                 "List of activation energies, $E_a$, for background material and compositional fields, "
                                 "for a total of N+1 values, where N is the number of compositional fields. "
                                 "If only one value is given, then all use the same value. "
                                 "Units: \\si{\\joule\\per\\mole}.");
              prm.declare_entry ("Activation volumes for diffusion creep", "6e-6",
                                 Patterns::Anything(),
                                 "List of activation volumes, $V_a$, for background material and compositional fields, "
                                 "for a total of N+1 values, where N is the number of compositional fields. "
                                 "If only one value is given, then all use the same value. "
                                 "Units: \\si{\\meter\\cubed\\per\\mole}.");
              prm.declare_entry ("Grain size", "1e-3", Patterns::Double (0.),
                                 "Units: \\si{\\meter}.");

              // parameters for dislocation creep
              prm.declare_entry ("Prefactors for dislocation creep", "1.1e-16",
                                 Patterns::Anything(),
                                 "List of viscosity prefactors, $A$, for background material and compositional fields, "
                                 "for a total of N+1 values, where N is the number of compositional fields. "
                                 "If only one value is given, then all use the same value. "
                                 "Units: \\si{\\pascal}$^{-n_{\\text{dislocation}}}$ \\si{\\per\\second}.");
              prm.declare_entry ("Stress exponents for dislocation creep", "3.5",
                                 Patterns::Anything(),
                                 "List of stress exponents, $n_{\\text{dislocation}}$, for background material and compositional fields, "
                                 "for a total of N+1 values, where N is the number of compositional fields. "
                                 "If only one value is given, then all use the same value.  Units: None.");
              prm.declare_entry ("Activation energies for dislocation creep", "530e3",
                                 Patterns::Anything(),
                                 "List of activation energies, $E_a$, for background material and compositional fields, "
                                 "for a total of N+1 values, where N is the number of compositional fields. "
                                 "If only one value is given, then all use the same value. "
                                 "Units: \\si{\\joule\\per\\mole}.");
              prm.declare_entry ("Activation volumes for dislocation creep", "1.4e-5",
                                 Patterns::Anything(),
                                 "List of activation volumes, $V_a$, for background material and compositional fields, "
                                 "for a total of N+1 values, where N is the number of compositional fields. "
                                 "If only one value is given, then all use the same value. "
                                 "Units: \\si{\\meter\\cubed\\per\\mole}.");

              // limits of viscosity
              prm.declare_entry ("Minimum viscosity", "1e17", Patterns::Double(0.),
                                 "Lower cutoff for effective viscosity. Units: \\si{\\pascal\\second}.");
              prm.declare_entry ("Maximum viscosity", "1e28", Patterns::Double(0.),
                                 "Upper cutoff for effective viscosity. Units: \\si{\\pascal\\second}.");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      DiffusionDislocation<dim>::parse_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Material model");
        {
          prm.enter_subsection("Viscosity");
          {
            prm.enter_subsection("Diffusion dislocation");
            {
              // parameters for diffusion creep
              const unsigned int n_fields = this->n_compositional_fields() + 1;

              A_diff = Utilities::possibly_extend_from_1_to_N(
                Utilities::string_to_double(Utilities::split_string_list(prm.get("Prefactors for diffusion creep"))),
                n_fields,
                "Prefactors for diffusion creep");

              m_diff = Utilities::possibly_extend_from_1_to_N(
                Utilities::string_to_double(Utilities::split_string_list(prm.get("Grain size exponents for diffusion creep"))),
  n_fields,
                "Grain size exponents for diffusion creep");

              E_diff = Utilities::possibly_extend_from_1_to_N(
                Utilities::string_to_double(Utilities::split_string_list(prm.get("Activation energies for diffusion creep"))),
                n_fields,
                "Activation energies for diffusion creep");

              V_diff = Utilities::possibly_extend_from_1_to_N(
                Utilities::string_to_double(Utilities::split_string_list(prm.get("Activation volumes for diffusion creep"))),
                n_fields,
                "Activation volumes for diffusion creep");

              d = prm.get_double("Grain size");

              // parameters for dislocation creep
              A_disl = Utilities::possibly_extend_from_1_to_N(
                Utilities::string_to_double(Utilities::split_string_list(prm.get("Prefactors for dislocation creep"))),
                n_fields,
                "Prefactors for dislocation creep");

              n_disl = Utilities::possibly_extend_from_1_to_N(
                Utilities::string_to_double(Utilities::split_string_list(prm.get("Stress exponents for dislocation creep"))),
                n_fields,
                "Stress exponents for dislocation creep");

              E_disl = Utilities::possibly_extend_from_1_to_N(
                Utilities::string_to_double(Utilities::split_string_list(prm.get("Activation energies for dislocation creep"))),
                n_fields,
                "Activation energies for dislocation creep");

              V_disl = Utilities::possibly_extend_from_1_to_N(
                Utilities::string_to_double(Utilities::split_string_list(prm.get("Activation volumes for dislocation creep"))),
                n_fields,
                "Activation volumes for dislocation creep");

              // limits of viscosity
              eta_min = prm.get_double("Minimum viscosity");
              eta_max = prm.get_double("Maximum viscosity");

              AssertThrow(eta_min > 0,
                          ExcMessage("The minimum viscosity must be greater than 0."));
              AssertThrow(eta_max > eta_min,
                          ExcMessage("The maximum viscosity must be greater than the minimum viscosity"));
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }

    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace MaterialModel
  {
    namespace Viscosity
    {
      ELASPECT_REGISTER_VISCOSITY_MODEL(DiffusionDislocation,
                                        "diffusion dislocation",
                                        "")
    }
  }
}
