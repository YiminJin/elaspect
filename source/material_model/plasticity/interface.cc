#include <elaspect/material_model/plasticity/interface.h>

#include <tuple>

namespace elaspect
{
  namespace MaterialModel
  {
    namespace Plasticity
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


// -------------------------------- Deal with registering plasticity models and automating
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
      register_plasticity_model (const std::string &name,
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
      create_plasticity_model (ParameterHandler &prm)
      {
        std::string model_name;
        prm.enter_subsection ("Material model");
        {
          prm.enter_subsection ("Plasticity");
          {
            model_name = prm.get ("Model name");
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();

        // We allow the plasticity model to be unspecified.
        Interface<dim> *plugin = nullptr;
        if (model_name != "unspecified")
          plugin = std::get<dim>(registered_plugins).create_plugin(
            model_name, "Material model::Plasticity::Model name");

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
          prm.enter_subsection ("Plasticity");
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


// explicit instantiations
namespace elaspect
{
  namespace internal
  {
    namespace Plugins
    {
      template <>
      std::list<internal::Plugins::PluginList<MaterialModel::Plasticity::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<MaterialModel::Plasticity::Interface<2> >::plugins = nullptr;

      template <>
      std::list<internal::Plugins::PluginList<MaterialModel::Plasticity::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<MaterialModel::Plasticity::Interface<3> >::plugins = nullptr;
    }
  }

  namespace MaterialModel
  {
    namespace Plasticity
    {
#define INSTANTIATE(dim) \
      template class Interface<dim>; \
      \
      template \
      void \
      register_plasticity_model<dim> (const std::string &, \
                                      const std::string &, \
                                      void ( *) (ParameterHandler &), \
                                      Interface<dim> *( *) ()); \
      \
      template \
      Interface<dim> * \
      create_plasticity_model<dim> (ParameterHandler &); \
      \
      template  \
      void \
      declare_parameters<dim> (ParameterHandler &);

      ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}
