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


#ifndef _elaspect_simulator_signals_h
#define _elaspect_simulator_signals_h

#include <elaspect/simulator_access.h>

#include <boost/signals2.hpp>

namespace elaspect
{
  namespace Assemblers
  {
    template <int dim>
    class Manager;
  }

  template <int dim>
  struct SimulatorSignals
  {
    /**
     * A signal that is called before setting up the initial conditions.
     *
     * The functions (slots) that can attach to this signal need to take one
     * argument: A SimulatorAccess object that describes the simulator.
     */
    boost::signals2::signal<void (parallel::distributed::Triangulation<dim> &)>  pre_set_initial_state;

    /**
     * A signal that is triggered at the end of the set_assemblers() function that
     * allows modification of the assembly objects active in this simulation.
     */
    boost::signals2::signal<void (const SimulatorAccess<dim> &,
                                  Assemblers::Manager<dim> &)> set_assemblers;
  };

  namespace internals
  {
    namespace SimulatorSignals
    {
      /**
       * Two functions that (in 2d and 3d) put a user-provided function onto a list
       * of functions that the Simulator object will later go through when
       * letting plugins connect their slots to signals.
       */
      void register_connector_function_2d(const std::function<void (elaspect::SimulatorSignals<2> &)> &connector);

      void register_connector_function_3d(const std::function<void (elaspect::SimulatorSignals<3> &)> &connector);

      /**
       * A function that is called by the Simulator object and that goes
       * through the list (with the corresponding dimension) created by the
       * previous pair of functions and call each of the user-provided
       * connector functions to let them register their slots with the
       * corresponding signals.
       */
      template <int dim>
      void call_connector_functions(elaspect::SimulatorSignals<dim> &signals);
    }
  }


  /**
   * A macro that is used in user-provided plugins to register a function that
   * is called at the beginning of a simulation by a Simulator object. When called,
   * the provided function will receive a SimulatorSignals object that contains
   * signals to which one can subscribe.
   *
   * For technical reasons, the macro takes two arguments denoting functions for the
   * 2d and 3d cases. These can, for example, be the names of 2d and 3d
   * instantiations of the same template function.
   */
#define ELASPECT_REGISTER_SIGNALS_CONNECTOR(connector_function_2d,connector_function_3d) \
  namespace ELASPECT_REGISTER_SIGNALS_CONNECTOR \
  { \
    struct dummy_do_register \
    { \
      dummy_do_register () \
      { \
        elaspect::internals::SimulatorSignals::register_connector_function_2d(connector_function_2d); \
        elaspect::internals::SimulatorSignals::register_connector_function_3d(connector_function_3d); \
      } \
    } dummy_variable; \
  }

  /**
   * A macro that is used to register a function that can be used to connect user
   * extension functions to the parameter-related signals declared in SimulatorSignals.
   *
   * In essence, this function simply registers a (global) function that is called
   * at the beginning of the program and that can be used to connect parameter
   * declaration and parsing functions to the signals listed above.
   */
#define ELASPECT_REGISTER_SIGNALS_PARAMETER_CONNECTOR(connector_function) \
  namespace ELASPECT_REGISTER_SIGNALS_PARAMETER_CONNECTOR_ ## connector_function \
  { \
    struct dummy_do_register_ ## connector_function \
    { \
      dummy_do_register_ ## connector_function () \
      { \
        connector_function (); \
      } \
    } dummy_variable_ ## classname; \
  }
}

#endif
