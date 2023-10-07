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


#ifndef _elaspect_termination_criteria_interface_h
#define _elaspect_termination_criteria_interface_h

#include <elaspect/plugins.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace TerminationCriteria
  {
    template <int dim>
    class Interface
    {
      public:
        virtual ~Interface ();

        virtual void initialize ();

        virtual bool execute () = 0;

        virtual double check_for_last_time_step (const double time_step) const;

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);
    };


    template <int dim>
    class Manager : public SimulatorAccess<dim>
    {
      public:
        bool execute () const;

        double check_for_last_time_step (const double time_step) const;

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);

        static
        void
        register_termination_criterion (const std::string &name,
                                        const std::string &description,
                                        void (*declare_parameters_function) (ParameterHandler &),
                                        Interface<dim> *(*factory_function) ());

      private:
        std::list<std::unique_ptr<Interface<dim> > > termination_objects;

        std::list<std::string>                       termination_obj_names;
    };


    /**
     * Given a class name, a name, and a description for the parameter file
     * for a termination criterion object, register it with the
     * elaspect::TerminationCriteria::Manager class.
     *
     * @ingroup TerminationCriteria
     */
#define ELASPECT_REGISTER_TERMINATION_CRITERION(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ELASPECT_REGISTER_TERMINATION_CRITERION_ ## classname \
  { \
    elaspect::internal::Plugins::RegisterHelper<elaspect::TerminationCriteria::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&elaspect::TerminationCriteria::Manager<2>::register_termination_criterion, \
                                name, description); \
    elaspect::internal::Plugins::RegisterHelper<elaspect::TerminationCriteria::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&elaspect::TerminationCriteria::Manager<3>::register_termination_criterion, \
                                name, description); \
  }
  }
}

#endif
