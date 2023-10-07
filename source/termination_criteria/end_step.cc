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


#include <elaspect/termination_criteria/end_step.h>

namespace elaspect
{
  namespace TerminationCriteria
  {
    template <int dim>
    bool
    EndStep<dim>::execute()
    {
      return (this->get_timestep_number () > end_step);
    }

    template <int dim>
    void
    EndStep<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Termination criteria");
      {
        prm.declare_entry ("End step", "100",
                           Patterns::Integer (0),
                           "Terminate the simulation once the specified timestep has been reached.");
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    EndStep<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Termination criteria");
      {
        end_step = prm.get_integer ("End step");
      }
      prm.leave_subsection ();
    }
  }
}

// explicit instantiations
namespace elaspect
{
  namespace TerminationCriteria
  {
    ELASPECT_REGISTER_TERMINATION_CRITERION(EndStep,
                                        "end step",
                                        "Terminate the simulation once the specified timestep "
                                        "has been reached. ")
  }
}
