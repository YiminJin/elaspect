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


#ifndef _elaspect_time_stepping_interface_h
#define _elaspect_time_stepping_interface_h

#include <elaspect/simulator_access.h>
#include <elaspect/termination_criteria/interface.h>

namespace elaspect
{
  namespace TimeStepping
  {
    using namespace dealii;

    template <int dim>
    class Interface
    {
      public:
        virtual ~Interface ();

        virtual void initialize ();

        virtual void update ();

        virtual double get_next_time_step_size () const = 0;

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
        void update ();

        double get_next_time_step_size () const;

        double get_minimum_time_step_size () const;

        double get_maximum_time_step_size () const;

        bool should_refine_mesh () const;

        bool should_simulation_terminate_now () const;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm);

        static
        void
        register_time_stepping_model (const std::string &name,
                                      const std::string &description,
                                      void (*declare_parameters_function) (ParameterHandler &),
                                      Interface<dim> *(*factory_function) ());

      private:
        double next_time_step_size;

        double first_time_step_size;

        double minimum_time_step_size;

        double maximum_time_step_size;

        TerminationCriteria::Manager<dim> termination_manager;

        std::list<std::unique_ptr<Interface<dim> > > plugin_list;
    };


#define ELASPECT_REGISTER_TIME_STEPPING_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ELASPECT_REGISTER_TIME_STEPPING_MODEL_ ## classname \
  { \
    elaspect::internal::Plugins::RegisterHelper<elaspect::TimeStepping::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&elaspect::TimeStepping::Manager<2>::register_time_stepping_model, \
                                name, description); \
    elaspect::internal::Plugins::RegisterHelper<elaspect::TimeStepping::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&elaspect::TimeStepping::Manager<3>::register_time_stepping_model, \
                                name, description); \
  }
  }
}

#endif
