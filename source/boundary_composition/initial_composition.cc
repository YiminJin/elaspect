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


#include <elaspect/boundary_composition/initial_composition.h>

namespace elaspect
{
  namespace BoundaryComposition
  {
    template <int dim>
    double
    InitialComposition<dim>::
    boundary_composition(const types::boundary_id /*boundary_indicator*/,
                         const Point<dim> &position,
                         const unsigned int compositional_field) const
    {
      return this->get_initial_composition_manager().initial_composition(position,
                                                                         compositional_field);
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace BoundaryComposition
  {
    ELASPECT_REGISTER_BOUNDARY_COMPOSITION_MODEL(InitialComposition,
                                             "initial composition",
                                             "")
  }
}
