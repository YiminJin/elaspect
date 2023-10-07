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


#ifndef _elaspect_boundary_composition_initial_composition_h
#define _elaspect_boundary_composition_initial_composition_h

#include <elaspect/initial_composition/interface.h>
#include <elaspect/boundary_composition/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace BoundaryComposition
  {
    template <int dim>
    class InitialComposition : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        double boundary_composition (const types::boundary_id boundary_indicator,
                                     const Point<dim> &position,
                                     const unsigned int compositional_field) const override;
    };
  }
}

#endif
