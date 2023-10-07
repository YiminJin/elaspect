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


#include <elaspect/time_stepping/inherited_time_step.h>

namespace elaspect
{
  namespace TimeStepping
  {
    template <int dim>
    double
    InheritedTimeStep<dim>::get_next_time_step_size() const
    {
      if (this->get_timestep_number() == 0)
        return std::numeric_limits<double>::max();
      else
        return this->get_timestep();
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace TimeStepping
  {
    ELASPECT_REGISTER_TIME_STEPPING_MODEL(InheritedTimeStep,
                                          "inherited time step",
                                          "")
  }
}
