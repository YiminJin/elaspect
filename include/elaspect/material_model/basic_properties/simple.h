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


#ifndef _elaspect_material_model_basic_properties_simple_h
#define _elaspect_material_model_basic_properties_simple_h

#include <elaspect/material_model/basic_properties/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace MaterialModel
  {
    namespace BasicProperties
    {
      using namespace dealii;

      template <int dim>
      class Simple : public Interface<dim>, public SimulatorAccess<dim>
      {
        public:
          void
          compute_basic_properties (const MaterialModel::MaterialModelInputs<dim> &in,
                                    const unsigned int i,
                                    MaterialModel::BasicPropertyOutputs<dim> &out) const override;

          FieldDependences::Dependence
          get_field_dependences (const MaterialProperties::Property requested_properties) const override;

          static
          void 
          declare_parameters (ParameterHandler &prm);

          void
          parse_parameters (ParameterHandler &prm) override;

        private:
          double reference_temperature;
          std::vector<double> reference_densities;
          std::vector<double> specific_heat_capacities;
          std::vector<double> thermal_expansivities;
          std::vector<double> thermal_conductivities;
          std::vector<double> bulk_moduli;
          std::vector<double> shear_moduli;
      };
    }
  }
}

#endif
