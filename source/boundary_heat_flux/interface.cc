#include <elaspect/boundary_heat_flux/interface.h>

namespace elaspect
{
  namespace BoundaryHeatFlux
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
    void
    Interface<dim>::
    declare_parameters (ParameterHandler &/*prm*/)
    {}


    template <int dim>
    void
    Interface<dim>::parse_parameters (ParameterHandler &/*prm*/)
    {}


// -------------------------------- Deal with registering models and automating
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
    register_boundary_heat_flux (const std::string &name,
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
    create_boundary_heat_flux (ParameterHandler &prm)
    {
      std::string model_name;
      prm.enter_subsection ("Boundary heat flux model");
      {
        model_name = prm.get ("Model name");
      }
      prm.leave_subsection ();

      return std::get<dim>(registered_plugins).create_plugin (model_name,
                                                              "Boundary heat flux model::Model name");
    }


    template <int dim>
    void
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Boundary heat flux model");
      const std::string pattern_of_names
        = std::get<dim>(registered_plugins).get_pattern_of_names ();
      prm.declare_entry ("Model name", "function",
                         Patterns::Selection (pattern_of_names),
                         "Select one of the following plugins:\n\n"
                         +
                         std::get<dim>(registered_plugins).get_description_string());
      prm.leave_subsection ();

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
      std::list<internal::Plugins::PluginList<BoundaryHeatFlux::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<BoundaryHeatFlux::Interface<2> >::plugins = nullptr;
      template <>
      std::list<internal::Plugins::PluginList<BoundaryHeatFlux::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<BoundaryHeatFlux::Interface<3> >::plugins = nullptr;
    }
  }

  namespace BoundaryHeatFlux
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  \
  template \
  void \
  register_boundary_heat_flux<dim> (const std::string &, \
                                    const std::string &, \
                                    void ( *) (ParameterHandler &), \
                                    Interface<dim> *( *) ()); \
  \
  template  \
  void \
  declare_parameters<dim> (ParameterHandler &); \
  \
  template \
  Interface<dim> * \
  create_boundary_heat_flux<dim> (ParameterHandler &prm);

    ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
