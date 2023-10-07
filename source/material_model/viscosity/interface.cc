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


#include <elaspect/material_model/viscosity/interface.h>

#include <deal.II/base/signaling_nan.h>

#include <tuple>

namespace elaspect
{
  namespace MaterialModel
  {
    namespace Viscosity
    {
      template <int dim>
      Interface<dim>::~Interface ()
      {}


      template <int dim>
      void
      Interface<dim>::initialize ()
      {}


      template <int dim>
      void
      Interface<dim>::update ()
      {}


      template <int dim>
      FieldDependences::Dependence
      Interface<dim>::get_field_dependences () const
      {
        return FieldDependences::none;
      }


      template <int dim>
      void
      Interface<dim>::declare_parameters (ParameterHandler &)
      {}


      template <int dim>
      void
      Interface<dim>::parse_parameters (ParameterHandler &)
      {}


// -------------------------------- Deal with registering viscosity models and automating
// -------------------------------- their setup and selection at run time

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
      register_viscosity_model (const std::string &name,
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
      Interface<dim> *
      create_viscosity_model (ParameterHandler &prm)
      {
        std::string model_name;
        prm.enter_subsection ("Material model");
        {
          prm.enter_subsection ("Viscosity");
          {
            model_name = prm.get ("Model name");
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();

        // We allow the viscosity model to be unspecified.
        Interface<dim> *plugin = nullptr;
        if (model_name != "unspecified")
          plugin = std::get<dim>(registered_plugins).create_plugin(
            model_name, "Material model::Viscosity::Model name");

        return plugin;
      }


      template <int dim>
      std::string
      get_valid_model_names_pattern ()
      {
        return std::get<dim>(registered_plugins).get_pattern_of_names ();
      }


      template <int dim>
      void
      declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection ("Material model");
        {
          prm.enter_subsection ("Viscosity");
          {
            const std::string pattern_of_names = get_valid_model_names_pattern<dim>();
            prm.declare_entry ("Model name", "unspecified",
                               Patterns::Selection (pattern_of_names+"|unspecified"),
                               "Select one of the following models:\n\n"
                               +
                               std::get<dim>(registered_plugins).get_description_string());
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();

        std::get<dim>(registered_plugins).declare_parameters (prm);
      }

    }
  }
}

namespace elaspect
{
  namespace internal
  {
    namespace Plugins
    {
      template <>
      std::list<internal::Plugins::PluginList<MaterialModel::Viscosity::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<MaterialModel::Viscosity::Interface<2> >::plugins = nullptr;

      template <>
      std::list<internal::Plugins::PluginList<MaterialModel::Viscosity::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<MaterialModel::Viscosity::Interface<3> >::plugins = nullptr;
    }
  }

  namespace MaterialModel
  {
    namespace Viscosity
    {
#define INSTANTIATE(dim) \
      template class Interface<dim>; \
      \
      template \
      void \
      register_viscosity_model<dim> (const std::string &, \
                                     const std::string &, \
                                     void ( *) (ParameterHandler &), \
                                     Interface<dim> *( *) ()); \
      \
      template \
      Interface<dim> * \
      create_viscosity_model<dim> (ParameterHandler &); \
      \
      template  \
      void \
      declare_parameters<dim> (ParameterHandler &);

    ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}
