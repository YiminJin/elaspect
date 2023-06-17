#include <elaspect/material_model/basic_properties/interface.h>

#include <deal.II/base/signaling_nan.h>

#include <tuple>

namespace elaspect
{
  namespace MaterialModel
  {
    template <int dim>
    BasicPropertyOutputs<dim>::BasicPropertyOutputs (const unsigned int n_fields)
      :
      specific_heat_capacities (n_fields, numbers::signaling_nan<double>()),
      thermal_expansivities (n_fields, numbers::signaling_nan<double>()),
      thermal_conductivities (n_fields, numbers::signaling_nan<double>()),
      densities (n_fields, numbers::signaling_nan<double>()),
      reaction_terms (n_fields, numbers::signaling_nan<double>()),
      bulk_moduli (n_fields, numbers::signaling_nan<double>()),
      shear_moduli (n_fields, numbers::signaling_nan<double>())
    {}


    namespace BasicProperties
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
      Interface<dim>::
      get_field_dependences (const MaterialProperties::Property /*requested_properties*/) const
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


// -------------------------------- Deal with registering basic property models and automating
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
      register_basic_property_model (const std::string &name,
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
      create_basic_property_model (ParameterHandler &prm)
      {
        std::string model_name;

        prm.enter_subsection ("Material model");
        {
          prm.enter_subsection ("Basic properties");
          {
            model_name = prm.get ("Model name");
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();

        // If one sets the model name to an empty string in the input file,
        // ParameterHandler produces an error while reading the file. However,
        // if one omits specifying any model name at all (not even setting it to
        // the empty string) then the value we get here is the empty string. If
        // we don't catch this case here, we end up with awkward downstream
        // errors because the value obviously does not conform to the Pattern.
        AssertThrow(model_name != "unspecified",
                    ExcMessage("You need to select a basic material model "
                               "(`set Model name' in `subsection Material model::Basic properties')."));

        Interface<dim> *plugin = 
          std::get<dim>(registered_plugins).create_plugin (model_name,
                                                           "Material model::Basic properties::Model name");

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
          prm.enter_subsection ("Basic properties");
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
      std::list<internal::Plugins::PluginList<MaterialModel::BasicProperties::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<MaterialModel::BasicProperties::Interface<2> >::plugins = nullptr;

      template <>
      std::list<internal::Plugins::PluginList<MaterialModel::BasicProperties::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<MaterialModel::BasicProperties::Interface<3> >::plugins = nullptr;
    }
  }

  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
    template struct BasicPropertyOutputs<dim>; \
    \
    namespace BasicProperties \
    { \
      template class Interface<dim>; \
      \
      template \
      void \
      register_basic_property_model<dim> (const std::string &, \
                                          const std::string &, \
                                          void ( *) (ParameterHandler &), \
                                          Interface<dim> *( *) ()); \
      \
      template \
      Interface<dim> * \
      create_basic_property_model<dim> (ParameterHandler &); \
      \
      template  \
      void \
      declare_parameters<dim> (ParameterHandler &); \
    }

    ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
