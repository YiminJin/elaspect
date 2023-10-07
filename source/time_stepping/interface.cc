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


#include <elaspect/time_stepping/interface.h>
#include <elaspect/simulator.h>

namespace elaspect
{
  namespace TimeStepping
  {
    /*----------------------------- Interface -------------------------------*/

    template <int dim>
    Interface<dim>::~Interface ()
    {}


    template <int dim>
    void Interface<dim>::initialize ()
    {}


    template <int dim>
    void Interface<dim>::update ()
    {}


    template <int dim>
    void 
    Interface<dim>::declare_parameters (ParameterHandler &)
    {}


    template <int dim>
    void 
    Interface<dim>::parse_parameters (ParameterHandler &)
    {}


    /*----------------------------- Manager -------------------------------*/

    template <int dim>
    void Manager<dim>::update ()
    {
      double new_time_step = std::numeric_limits<double>::max();
      for (const auto &p : plugin_list)
      {
        p->update();
        new_time_step = std::min(new_time_step, p->get_next_time_step_size());
      }

      if (this->get_timestep_number() == 0)
        new_time_step = std::min(first_time_step_size, new_time_step);

      new_time_step = Utilities::MPI::min(new_time_step, this->get_mpi_communicator());

      // Do not go below a minimum time step size as given by the user. Note that
      // we first apply the minimum and then apply further restrictions that might
      // lower this value. Hitting the end time is more important than a minimum
      // time step size, for example.
      new_time_step = std::max(minimum_time_step_size, new_time_step);

      // Make sure we do not exceed the maximum time step length.
      new_time_step = std::min(new_time_step, maximum_time_step_size);

      // Make sure we reduce the time step length appropriately if we terminate after this step
      next_time_step_size = std::min(new_time_step, termination_manager.check_for_last_time_step(new_time_step));
    }


    template <int dim>
    double
    Manager<dim>::get_next_time_step_size () const
    {
      return next_time_step_size;
    }


    template <int dim>
    bool
    Manager<dim>::should_simulation_terminate_now () const
    {
      return termination_manager.execute();
    }


    template <int dim>
    double Manager<dim>::get_maximum_time_step_size() const
    {
      return maximum_time_step_size;
    }


    template <int dim>
    double Manager<dim>::get_minimum_time_step_size() const
    {
      return minimum_time_step_size;
    }


    namespace
    {
      std::tuple
      <void *,
      void *,
      elaspect::internal::Plugins::PluginList<Interface<2> >,
      elaspect::internal::Plugins::PluginList<Interface<3> > > registered_plugins;
    }


    template <int dim>
    void
    Manager<dim>::register_time_stepping_model(const std::string &name,
                                               const std::string &description,
                                               void (*declare_parameters_function) (ParameterHandler &),
                                               Interface<dim> *(*factory_function) ())
    {
      std::get<dim>(registered_plugins).register_plugin (name,
                                                         description,
                                                         declare_parameters_function,
                                                         factory_function);
    }


    template <int dim>
    void
    Manager<dim>::declare_parameters (ParameterHandler &prm)
    {
      TerminationCriteria::Manager<dim>::declare_parameters(prm);

      prm.enter_subsection("Time stepping");
      {
        prm.declare_entry("First time step size",
                          /* boost::lexical_cast<std::string>(std::numeric_limits<double>::max()) /
                                                              year in seconds = */ "5.69e+300",
                          Patterns::Double(0),
                          "Specify the first time step size, since the velocity field is unknown "
                          "at the first time step.");
        prm.declare_entry("Minimum time step size",
                          /* boost::lexical_cast<std::string>(std::numeric_limits<double>::min())
                             = */ "2.23e-308",
                          Patterns::Double(0),
                          "Specify a minimum time step size.");

        prm.declare_entry("Maximum time step size", 
                          /* boost::lexical_cast<std::string>(std::numeric_limits<double>::max()) /
                                                              year in seconds = */ "5.69e+300",
                          Patterns::Double(0),
                          "Specify a maximum time step size.");

        const std::string pattern_of_names
          = std::get<dim>(registered_plugins).get_pattern_of_names ();

        prm.declare_entry("List of model names", "",
                          Patterns::MultipleSelection(pattern_of_names),
                          "A comma separated list of time stepping plugins that "
                          "will be used to calculate the time step size. The minimum of the "
                          "result of each plugin will be used.\n\n"
                          "The following plugins are available:\n\n"
                          +
                          std::get<dim>(registered_plugins).get_description_string());
      }
      prm.leave_subsection();

      std::get<dim>(registered_plugins).declare_parameters(prm);
    }


    template <int dim>
    void
    Manager<dim>::parse_parameters (ParameterHandler &prm)
    {
      termination_manager.initialize_simulator(this->get_simulator());
      termination_manager.parse_parameters(prm);

      std::vector<std::string> model_names;
      prm.enter_subsection("Time stepping");
      {
        first_time_step_size   = prm.get_double("First time step size")
                                 * (this->convert_output_to_years() ? year_in_seconds : 1.0);
        minimum_time_step_size = prm.get_double("Minimum time step size")
                                 * (this->convert_output_to_years() ? year_in_seconds : 1.0);
        maximum_time_step_size = prm.get_double("Maximum time step size")
                                 * (this->convert_output_to_years() ? year_in_seconds : 1.0);

        model_names = Utilities::split_string_list(prm.get("List of model names"));

        if (model_names.size() == 0)
          model_names.emplace_back("convection time step");

        AssertThrow(Utilities::has_unique_entries(model_names),
                    ExcMessage("The list of strings for the parameter "
                               "'Time stepping/List of model names' contains entries more than once. "
                               "This is not allowed. Please check your parameter file."));

        if (this->get_parameters().use_ALE_method)
        {
          AssertThrow(std::find(model_names.begin(), model_names.end(), "convection time step")
                      != model_names.end(),
                      ExcMessage("When ALE method is applied, 'convection time step' must be "
                                 "included in 'Time stepping/List of model names'."));
        }
      }
      prm.leave_subsection();

      for (const auto &plugin_name : model_names)
      {
        plugin_list.push_back (std::unique_ptr<Interface<dim> >
                               (std::get<dim>(registered_plugins)
                                .create_plugin (plugin_name, 
                                                "Time stepping::Model names")));

        if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(&*plugin_list.back()))
          sim->initialize_simulator (this->get_simulator());

        plugin_list.back()->parse_parameters (prm);
        plugin_list.back()->initialize ();
      }
    }
  }
}


// explicit instantiations
namespace elaspect
{

  namespace internal
  {
    namespace Plugins
    {
      template <>
      std::list<internal::Plugins::PluginList<TimeStepping::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<TimeStepping::Interface<2> >::plugins = nullptr;
      template <>
      std::list<internal::Plugins::PluginList<TimeStepping::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<TimeStepping::Interface<3> >::plugins = nullptr;
    }
  }

  namespace TimeStepping
  {
#define INSTANTIATE(dim) \
    \
    template class Interface<dim>; \
    template class Manager<dim>; \

    ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
