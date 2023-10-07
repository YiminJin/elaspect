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


#include <elaspect/boundary_temperature/initial_temperature.h>
#include <elaspect/initial_temperature/interface.h>

namespace elaspect
{
  namespace BoundaryTemperature
  {
    template <int dim>
    double
    InitialTemperature<dim>::
    boundary_temperature(const types::boundary_id /*boundary_indicator*/,
                         const Point<dim> &position) const
    {
      return this->get_initial_temperature_manager().initial_temperature(position);
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace BoundaryTemperature
  {
    ELASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL(InitialTemperature,
                                                 "initial temperature",
                                                 "")
  }
}
