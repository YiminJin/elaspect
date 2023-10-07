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


#include <elaspect/simulator_signals.h>

namespace elaspect
{
  namespace internals
  {
    namespace SimulatorSignals
    {
      std::list<std::function<void (elaspect::SimulatorSignals<2> &)>> connector_functions_2d;
      std::list<std::function<void (elaspect::SimulatorSignals<3> &)>> connector_functions_3d;

      static bool connector_functions_have_been_called = false;

      // add a user-provided connector to the list of connectors we keep
      void register_connector_function_2d(const std::function<void (elaspect::SimulatorSignals<2> &)> &connector)
      {
        Assert(!connector_functions_have_been_called,
               ExcMessage("Registration of signal connector happened after connection has already been called!"));
        connector_functions_2d.push_back (connector);       
      }


      void register_connector_function_3d(const std::function<void (elaspect::SimulatorSignals<3> &)> &connector)
      {
        Assert(!connector_functions_have_been_called,
               ExcMessage("Registration of signal connector happened after connection has already been called!"));
        connector_functions_3d.push_back (connector);
      }


      // call connectors to ensure that plugins get a change to register their slots
      template <>
      void call_connector_functions(elaspect::SimulatorSignals<2> &signals)
      {
        Assert(!connector_functions_have_been_called, 
               ExcInternalError());
        
        for (const auto &p : connector_functions_2d)
          p(signals);

        connector_functions_have_been_called = true;
      }


      template <>
      void call_connector_functions(elaspect::SimulatorSignals<3> &signals)
      {
        Assert(!connector_functions_have_been_called,
               ExcInternalError());

        for (const auto &p : connector_functions_3d)
          p(signals);

        connector_functions_have_been_called = true;
      }
    }
  }
}


// explicit instantiations
namespace elaspect
{
#define INSTANTIATE(dim) \
  template struct SimulatorSignals<dim>;


  ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
