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


#include <elaspect/postprocess/visualization/strain_rate.h>

namespace elaspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      StrainRate<dim>::StrainRate()
        :
        DataPostprocessorScalar<dim>("strain_rate", update_gradients)
      {}


      template <int dim>
      void
      StrainRate<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points, ExcInternalError());
        Assert (computed_quantities[0].size() == 1, ExcInternalError());
        Assert (input_data.solution_gradients[0].size() == this->introspection().n_components, ExcInternalError());

        for (unsigned int q = 0; q < n_quadrature_points; ++q)
        {
          Tensor<2,dim> grad_du;
          for (unsigned int d = 0; d < dim; ++d)
            grad_du[d] = input_data.solution_gradients[q][d];

          const SymmetricTensor<2,dim> depsilon = symmetrize(grad_du);
          computed_quantities[q](0) = std::sqrt(std::fabs(second_invariant(depsilon))) / this->get_timestep();
        }
      }
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ELASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(StrainRate,
                                                "strain rate",
                                                "")
    }
  }
}
