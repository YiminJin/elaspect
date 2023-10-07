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


#ifndef _elaspect_termination_criteria_end_time_h
#define _elaspect_termination_criteria_end_time_h

#include <elaspect/termination_criteria/interface.h>

namespace elaspect
{
  namespace TerminationCriteria
  {
    template <int dim>
    class EndTime : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        double 
        check_for_last_time_step (const double time_step) const override;
        
        bool
        execute () override;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        double end_time;
    };
  }
}

#endif
