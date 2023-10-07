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


#ifndef _elaspect_material_model_viscosity_simple_h
#define _elaspect_material_model_viscosity_simple_h

#include <elaspect/material_model/viscosity/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace MaterialModel
  {
    namespace Viscosity
    {
      using namespace dealii;

      template <int dim>
      class Simple : public Interface<dim>, public SimulatorAccess<dim>
      {
        public:
          double
          compute_viscosity (const MaterialModel::MaterialModelInputs<dim> &in,
                             const unsigned int i,
                             const unsigned int field) const override;

          FieldDependences::Dependence
          get_field_dependences () const override;

          static
          void
          declare_parameters (ParameterHandler &prm);

          virtual
          void
          parse_parameters (ParameterHandler &prm);

        private:
          double reference_T;
          std::vector<double> eta;
          std::vector<double> beta;
          double maximum_thermal_prefactor;
          double minimum_thermal_prefactor;
      };
    }
  }
}

#endif
