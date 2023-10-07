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


#include <elaspect/time_stepping/function.h>

namespace elaspect
{
  namespace TimeStepping
  {
    template <int dim>
    void Function<dim>::update ()
    {
      if (this->convert_output_to_years())
        time_stepping_function.set_time (this->get_time() / year_in_seconds);
      else
        time_stepping_function.set_time (this->get_time());
    }


    template <int dim>
    double
    Function<dim>::get_next_time_step_size () const
    {
      const double time_step_size = time_stepping_function.value(Point<dim>());
      if (this->convert_output_to_years())
        return time_step_size * year_in_seconds;
      else
        return time_step_size;
    }


    template <int dim>
    void
    Function<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Time stepping");
      {
        prm.enter_subsection("Function");
        {
          Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Function<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Time stepping");
      {
        prm.enter_subsection("Function");
        {
          try
          {
            time_stepping_function.parse_parameters(prm);
          }
          catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Time stepping.Function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'"
                      << "More information about the cause of the parse error \n"
                      << "is shown below.\n";
            throw;
          }
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace TimeStepping
  {
    ELASPECT_REGISTER_TIME_STEPPING_MODEL(Function,
                                      "function",
                                      "")
  }
}
