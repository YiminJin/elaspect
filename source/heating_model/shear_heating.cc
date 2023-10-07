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


#include <elaspect/heating_model/shear_heating.h>

namespace elaspect
{
  namespace HeatingModel
  {
    template <int dim>
    void
    ShearHeating<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
              const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      const unsigned int n_evaluation_points = material_model_inputs.n_evaluation_points();
      AssertDimension(heating_model_outputs.heating_source_terms.size(), n_evaluation_points);

      for (unsigned int q = 0; q < n_evaluation_points; ++q)
      {
        if (this->constitutive_relation() & ConstitutiveRelation::viscosity)
        {
          const SymmetricTensor<2,dim> deviatoric_strain_rate = 
            deviator(symmetrize(material_model_inputs.incremental_displacement_gradient[q]))
            / this->get_timestep();

          const SymmetricTensor<2,dim> viscous_stress =
            2.0 * material_model_outputs.viscosities[q] * deviatoric_strain_rate;

          heating_model_outputs.heating_source_terms[q] = viscous_stress * deviatoric_strain_rate;
        }
        else
          heating_model_outputs.heating_source_terms[q] = 0.0;
      }
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace HeatingModel
  {
    ELASPECT_REGISTER_HEATING_MODEL(ShearHeating,
                                "shear heating",
                                "Implementation of a standard model for shear heating. "
                                "Adds the term: "
                                "$  2 \\eta \\left( \\varepsilon - \\frac{1}{3} \\text{tr} "
                                "\\varepsilon \\mathbf 1 \\right) : \\left( \\varepsilon - \\frac{1}{3} "
                                "\\text{tr} \\varepsilon \\mathbf 1 \\right)$ to the "
                                "right-hand side of the temperature equation."
                                )
  }
}
