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


#include <elaspect/initial_topography/zero_topography.h>

namespace elaspect
{
  namespace InitialTopography
  {
    template <int dim>
    double
    ZeroTopography<dim>::
    value (const Point<dim-1> &/*p*/) const
    {
      // return a zero value regardless of position
      return 0.0;
    }


    template <int dim>
    double
    ZeroTopography<dim>::
    max_topography () const
    {
      return 0;
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace InitialTopography
  {
    ELASPECT_REGISTER_INITIAL_TOPOGRAPHY_MODEL(ZeroTopography,
                                           "zero topography",
                                           "Implementation of a model in which the initial topography "
                                           "is zero. ")
  }
}
