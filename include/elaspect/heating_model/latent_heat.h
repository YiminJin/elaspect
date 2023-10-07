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


#ifndef _elaspect_heating_model_latent_heat_h
#define _elaspect_heating_model_latent_heat_h

#include <elaspect/heating_model/interface.h>

namespace elaspect
{
  namespace HeatingModel
  {
    using namespace dealii;

    template <int dim>
    class LatentHeat : public Interface<dim>
    {
      public:
        /**
         * Compute the heating model outputs for this class.
         */
        void
        evaluate(const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                 const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                 HeatingModel::HeatingModelOutputs &heating_model_outputs) const override;
    };
  }
}

#endif
