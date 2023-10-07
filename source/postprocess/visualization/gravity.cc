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


#include <elaspect/postprocess/visualization/gravity.h>
#include <elaspect/gravity_model/interface.h>

namespace elaspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      Gravity<dim>::Gravity()
        :
        DataPostprocessorVector<dim>("gravity", update_quadrature_points)
      {}


      template <int dim>
      void
      Gravity<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.evaluation_points.size();
        Assert (computed_quantities.size() == n_quadrature_points, ExcInternalError());
        Assert (computed_quantities[0].size() == dim, ExcInternalError());

        for (unsigned int q = 0; q < n_quadrature_points; ++q)
        {
          const Tensor<1,dim> g = this->get_gravity_model().gravity_vector(input_data.evaluation_points[q]);
          for (unsigned int k = 0; k < dim; ++k)
            computed_quantities[q](k) = g[k];
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
      ELASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(Gravity,
                                                    "gravity",
                                                    "A visualization output object that outputs the gravity vector."
                                                    "\n\n"
                                                    "Physical units: \\si {\\meter\\per\\second\\squared} .")
    }
  }
}
