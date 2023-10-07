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


#include <elaspect/boundary_velocity/zero_velocity.h>


namespace elaspect
{
  namespace BoundaryVelocity
  {
    template <int dim>
    Tensor<1,dim>
    ZeroVelocity<dim>::
    boundary_velocity (const types::boundary_id,
                       const Point<dim> &) const
    {
      // return a zero tensor regardless of position
      return Tensor<1,dim>();
    }
  }
}

// explicit instantiations
namespace elaspect
{
  namespace BoundaryVelocity
  {
    ELASPECT_REGISTER_BOUNDARY_VELOCITY_MODEL(ZeroVelocity,
                                          "zero velocity",
                                          "Implementation of a model in which the boundary "
                                          "velocity is zero. This is commonly referred to as "
                                          "a ``stick boundary condition'', indicating that "
                                          "the material ``sticks'' to the material on the "
                                          "other side of the boundary.")
  }
}
