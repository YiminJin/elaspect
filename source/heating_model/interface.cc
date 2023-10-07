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


#include <elaspect/global.h>
#include <elaspect/utilities.h>
#include <elaspect/heating_model/interface.h>

namespace elaspect
{
  namespace HeatingModel
  {
    /*---------------------- HeatingModelOutputs ----------------------*/

    HeatingModelOutputs::HeatingModelOutputs(const unsigned int n_points)
      : heating_source_terms(n_points, numbers::signaling_nan<double>())
      , lhs_latent_heat_terms(n_points, numbers::signaling_nan<double>())
    {}


    void HeatingModelOutputs::reset()
    {
      for (unsigned int q = 0; q < heating_source_terms.size(); ++q)
      {
        heating_source_terms[q] = numbers::signaling_nan<double>();
        lhs_latent_heat_terms[q] = numbers::signaling_nan<double>();
      }
    }



    /*-------------------------- Interface --------------------------*/

    template <int dim>
    void
    Interface<dim>::initialize()
    {}


    template <int dim>
    void
    Interface<dim>::update()
    {}


    template <int dim>
    void
    Interface<dim>::declare_parameters(ParameterHandler &)
    {}


    template <int dim>
    void
    Interface<dim>::parse_parameters(ParameterHandler &)
    {}


    /*-------------------------- Manager --------------------------*/

    template <int dim>
    void Manager<dim>::update()
    {
      for (const auto &heating_model : heating_model_objects)
        heating_model->update();
    }


    template <int dim>
    void
    Manager<dim>::evaluate(const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                           const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                           HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      const unsigned int n_evaluation_points = material_model_inputs.n_evaluation_points();
      for (unsigned int q = 0; q < n_evaluation_points; ++q)
      {
        heating_model_outputs.heating_source_terms[q] = 0.0;
        heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
      }

      HeatingModel::HeatingModelOutputs individual_heating_outputs(n_evaluation_points);

      for (const auto &heating_model : heating_model_objects)
      {
        heating_model->evaluate(material_model_inputs, 
                                material_model_outputs, 
                                individual_heating_outputs);

        for (unsigned int q = 0; q < n_evaluation_points; ++q)
        {
          heating_model_outputs.heating_source_terms[q] += individual_heating_outputs.heating_source_terms[q];
          heating_model_outputs.lhs_latent_heat_terms[q] += individual_heating_outputs.lhs_latent_heat_terms[q];
          individual_heating_outputs.reset();
        }
      }
    }


    template <int dim>
    const std::vector<std::string> &
    Manager<dim>::get_active_heating_model_names() const
    {
      return model_names;
    }


    template <int dim>
    const std::list<std::unique_ptr<Interface<dim>>> &
    Manager<dim>::get_active_heating_models() const
    {
      return heating_model_objects;
    }


    namespace
    {
      std::tuple
      <void *,
      void *,
      elaspect::internal::Plugins::PluginList<Interface<2>>,
      elaspect::internal::Plugins::PluginList<Interface<3>>> registered_plugins;
    }


    template <int dim>
    void
    Manager<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Heating model");
      {
        const std::string pattern_of_names
          = std::get<dim>(registered_plugins).get_pattern_of_names ();

        prm.declare_entry("List of model names",
                          "",
                          Patterns::MultipleSelection(pattern_of_names),
                          "A comma separated list of heating models that "
                          "will be used to calculate the heating terms in the energy "
                          "equation. The results of each of these criteria, i.e., "
                          "the heating source terms and the latent heat terms for the "
                          "left hand side will be added.\n\n"
                          "The following heating models are available:\n\n"
                          +
                          std::get<dim>(registered_plugins).get_description_string());
      }
      prm.leave_subsection ();

      std::get<dim>(registered_plugins).declare_parameters (prm);
    }


    template <int dim>
    void
    Manager<dim>::register_heating_model (const std::string &name,
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
    Manager<dim>::parse_parameters(ParameterHandler &prm)
    {
      // find out which plugins are requested and the various other
      // parameters we declare here
      prm.enter_subsection("Heating model");
      {
        model_names = Utilities::split_string_list(prm.get("List of model names"));

        AssertThrow(Utilities::has_unique_entries(model_names),
                    ExcMessage("The list of strings for the parameter "
                               "'Heating model/List of model names' contains entries more than once. "
                               "This is not allowed. Please check your parameter file."));
      }
      prm.leave_subsection();

      // go through the list, create objects and let them parse
      // their own parameters
      for (auto &model_name : model_names)
      {
        heating_model_objects.push_back(std::unique_ptr<Interface<dim>>
                                        (std::get<dim>(registered_plugins)
                                         .create_plugin (model_name,
                                                         "Heating model::Model names")));

        if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(&*heating_model_objects.back()))
          sim->initialize_simulator (this->get_simulator());

        heating_model_objects.back()->parse_parameters (prm);
        heating_model_objects.back()->initialize ();
      }
    }


    template <int dim>
    std::string
    get_valid_model_names_pattern()
    {
      return std::get<dim>(registered_plugins).get_pattern_of_names();
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
      std::list<internal::Plugins::PluginList<HeatingModel::Interface<2>>::PluginInfo> *
      internal::Plugins::PluginList<HeatingModel::Interface<2>>::plugins = nullptr;
      template <>
      std::list<internal::Plugins::PluginList<HeatingModel::Interface<3>>::PluginInfo> *
      internal::Plugins::PluginList<HeatingModel::Interface<3>>::plugins = nullptr;
    }
  }

  namespace HeatingModel
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  template class Manager<dim>; \
  \
  template \
  std::string \
  get_valid_model_names_pattern<dim> ();

    ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
