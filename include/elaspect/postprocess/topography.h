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


#ifndef _elaspect_postprocess_topography_h
#define _elaspect_postprocess_topography_h

#include <elaspect/postprocess/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace Postprocess
  {
    template <int dim>
    class Topography : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        std::pair<std::string,std::string> execute (TableHandler &statistics) override;

        static
        void
        declare_parameters(ParameterHandler &prm);

        void
        parse_parameters(ParameterHandler &prm) override;

      private:
        bool write_to_file;

        double output_interval;

        double last_output_time;
    };
  }
}

#endif
