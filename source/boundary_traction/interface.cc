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


#include <elaspect/boundary_traction/interface.h>

namespace elaspect
{
  namespace BoundaryTraction
  {
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
    Interface<dim>::declare_parameters (dealii::ParameterHandler &)
    {}


    template <int dim>
    void
    Interface<dim>::parse_parameters (dealii::ParameterHandler &)
    {}


// -------------------------------- Deal with registering boundary_traction models and automating
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
    register_boundary_traction (const std::string &name,
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
    create_boundary_traction (const std::string &name)
    {
      Interface<dim> *plugin = std::get<dim>(registered_plugins).create_plugin (name,
                                                                                "Boundary traction conditions");
      return plugin;
    }


    template <int dim>
    std::string
    get_names ()
    {
      return std::get<dim>(registered_plugins).get_pattern_of_names ();
    }


    template <int dim>
    void
    declare_parameters (ParameterHandler &prm)
    {
      std::get<dim>(registered_plugins).declare_parameters (prm);
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
      std::list<internal::Plugins::PluginList<BoundaryTraction::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<BoundaryTraction::Interface<2> >::plugins = nullptr;
      template <>
      std::list<internal::Plugins::PluginList<BoundaryTraction::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<BoundaryTraction::Interface<3> >::plugins = nullptr;
    }
  }

  namespace BoundaryTraction
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  \
  template \
  void \
  register_boundary_traction<dim> (const std::string &, \
                                   const std::string &, \
                                   void ( *) (ParameterHandler &), \
                                   Interface<dim> *( *) ()); \
  \
  template  \
  void \
  declare_parameters<dim> (ParameterHandler &); \
  \
  template  \
  std::string \
  get_names<dim> (); \
  \
  template \
  Interface<dim> * \
  create_boundary_traction<dim> (const std::string &);

    ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
