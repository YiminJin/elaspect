/*
  Copyright (C) 2023 by Yimin Jin.

  This file is part of elASPECT.

  elASPECT is modified from the free software ASPECT; you can 
  redistribute it and/or modify it under the terms of the GNU 
  General Public License as published by the Free Software 
  Foundation; either version 2, or (at your option) any later 
  version.

  elASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with elASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <elaspect/material_model/viscosity/simple.h>
#include <elaspect/utilities.h>

namespace elaspect
{
  namespace MaterialModel
  {
    namespace Viscosity
    {
      template <int dim>
      double
      Simple<dim>::
      compute_viscosity(const MaterialModel::MaterialModelInputs<dim> &in,
                        const unsigned int i,
                        const unsigned int field) const
      {
        double temperature_dependence = 1.0;
        if (reference_T > 0)
        {
          const double dT = in.temperature[i] - reference_T;
          temperature_dependence = std::max(std::min(std::exp(-beta[field] * dT / reference_T),
                                                     maximum_thermal_prefactor),
                                            minimum_thermal_prefactor);
        }

        return temperature_dependence * eta[field];
      }


      template <int dim>
      FieldDependences::Dependence
      Simple<dim>::get_field_dependences() const
      {
        FieldDependences::Dependence dependence = FieldDependences::none;
        if (reference_T > 0)
          dependence |= FieldDependences::temperature;

        return dependence;
      }


      template <int dim>
      void
      Simple<dim>::declare_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Material model");
        {
          prm.enter_subsection("Viscosity");
          {
            prm.enter_subsection("Simple model");
            {
              prm.declare_entry("Reference temperature", "293.",
                                Patterns::Double(0.),
                                "The reference temperature $T_0$. Units: \\si{\\kelvin}.");
              prm.declare_entry("Viscosity", "5e24",
                                Patterns::List(Patterns::Double(0.)),
                                "List of of the constant viscosity, $\\eta_0$, for background "
                                "material and compositional fields, for a total of N+1 values, "
                                "where N is the number of compositional fields. If only one value "
                                "is given, then all use the same value. "
                                " Units: \\si{\\pascal\\second}.");
              prm.declare_entry("Thermal viscosity exponent", "0.0",
                                Patterns::List(Patterns::Double(0.)),
                                "List of the temperature dependence, $\\beta$, for background "
                                "material and compositional fields, for a total of N+1 values, "
                                "where N is the number of compositional fields. If only one value "
                                "is given, then all use the same value. ");
              prm.declare_entry("Maximum thermal prefactor", "1.0e2",
                                Patterns::Double(0.),
                                "The maximum value of the viscosity prefactor associated with "
                                "temperature dependence.");
              prm.declare_entry("Minimum thermal prefactor", "1.0e-2",
                                Patterns::Double(0.),
                                "The minimum value of the viscosity prefactor associated with "
                                "temperature dependence.");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      Simple<dim>::parse_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Material model");
        {
          prm.enter_subsection("Viscosity");
          {
            prm.enter_subsection("Simple model");
            {
              reference_T = prm.get_double("Reference temperature");

              const unsigned int n_fields = this->n_compositional_fields() + 1;

              eta = Utilities::possibly_extend_from_1_to_N(
                Utilities::string_to_double(Utilities::split_string_list(prm.get("Viscosity"))),
                n_fields,
                "Viscosity");

              beta = Utilities::possibly_extend_from_1_to_N(
                Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal viscosity exponent"))),
                n_fields,
                "Thermal viscosity exponent");

              maximum_thermal_prefactor = prm.get_double("Maximum thermal prefactor");
              minimum_thermal_prefactor = prm.get_double("Minimum thermal prefactor");

              if (maximum_thermal_prefactor == 0.0)
                maximum_thermal_prefactor = std::numeric_limits<double>::max();
              if (minimum_thermal_prefactor == 0.0)
                minimum_thermal_prefactor = std::numeric_limits<double>::min();

              if (reference_T == 0.0)
              {
                for (unsigned int i = 0; i < n_fields; ++i)
                  AssertThrow(beta[i] == 0, 
                              ExcMessage("Error: The simple viscosity model cannot take zero "
                                         "as reference temperature when temperature-dependent "
                                         "viscosity is taken into account."));
              }
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
      ELASPECT_REGISTER_VISCOSITY_MODEL(Simple,
                                    "simple",
                                    "")
    }
  }
}
