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


#include <elaspect/mesh_refinement/strain_rate.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace elaspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    StrainRate<dim>::execute (Vector<float> &indicators) const
    {
      indicators = 0;

      const QMidpoint<dim> quadrature;

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature,
                               update_quadrature_points | update_values | update_gradients);

      std::vector<SymmetricTensor<2,dim> > strain_increments (quadrature.size());

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
        {
          const unsigned int idx = cell->active_cell_index();
          fe_values.reinit(cell);

          fe_values[this->introspection().extractors.displacement].get_function_symmetric_gradients (
            this->get_solution(), strain_increments);

          indicators(idx) = strain_increments[0].norm();
        }
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace MeshRefinement
  {
    ELASPECT_REGISTER_MESH_REFINEMENT_CRITERION(StrainRate,
                                            "strain rate",
                                            "A mesh refinement criterion that computes the "
                                            "refinement indicators equal to the strain rate "
                                            "norm computed at the center of the elements.")
  }
}
