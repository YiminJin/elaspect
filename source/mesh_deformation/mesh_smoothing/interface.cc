#include <elaspect/mesh_deformation/mesh_smoothing/interface.h>

namespace elaspect
{
  namespace MeshDeformation
  {
    namespace MeshSmoothing
    {
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


// -------------------------------- Deal with registering mesh smoothing models and automating
// -------------------------------- their setup and selection at run time

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
      register_mesh_smoothing_model(const std::string &name,
                                    const std::string &description,
                                    void (*declare_parameters_function) (ParameterHandler &),
                                    Interface<dim> *(*factory_function) ())
      {
        std::get<dim>(registered_plugins).register_plugin(name,
                                                          description,
                                                          declare_parameters_function,
                                                          factory_function);
      }


      template <int dim>
      Interface<dim> *
      create_mesh_smoothing_model(ParameterHandler &prm)
      {
        std::string model_name;
        prm.enter_subsection("Mesh deformation");
        {
          prm.enter_subsection("Mesh smoothing");
          {
            model_name = prm.get("Model name");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();

        // If one sets the model name to an empty string in the input file,
        // ParameterHandler produces an error while reading the file. However,
        // if one omits specifying any model name at all (not even setting it to
        // the empty string) then the value we get here is the empty string. If
        // we don't catch this case here, we end up with awkward downstream
        // errors because the value obviously does not conform to the Pattern.
        AssertThrow(model_name != "unspecified",
                    ExcMessage("You need to select a mesh smoothing model "
                               "(`set Model name' in `subsection Mesh deformation::Mesh smoothing')."));
        return std::get<dim>(registered_plugins).create_plugin(model_name, 
                                                               "Mesh deformation::Mesh smoothing::Model name");
      }


      template <int dim>
      void
      declare_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Mesh deformation");
{
          prm.enter_subsection("Mesh smoothing");
          {
            const std::string pattern_of_names =
              std::get<dim>(registered_plugins).get_pattern_of_names();

            prm.declare_entry("Model name", "laplacian",
                              Patterns::Selection(pattern_of_names + "|unspecified"),
                              "Select one of the following models:\n\n"
                              +
                              std::get<dim>(registered_plugins).get_description_string());
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();

        std::get<dim>(registered_plugins).declare_parameters(prm);
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
      std::list<internal::Plugins::PluginList<MeshDeformation::MeshSmoothing::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<MeshDeformation::MeshSmoothing::Interface<2> >::plugins = nullptr;

      template <>
      std::list<internal::Plugins::PluginList<MeshDeformation::MeshSmoothing::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<MeshDeformation::MeshSmoothing::Interface<3> >::plugins = nullptr;
    }
  }

  namespace MeshDeformation
  {
    namespace MeshSmoothing
    {
#define INSTANTIATE(dim) \
      template class Interface<dim>; \
      \
      template \
      void \
      register_mesh_smoothing_model<dim>(const std::string &, \
                                         const std::string &, \
                                         void ( *) (ParameterHandler &), \
                                         Interface<dim> *( *) ()); \
      \
      template \
      Interface<dim> * \
      create_mesh_smoothing_model<dim>(ParameterHandler &); \
      \
      template  \
      void \
      declare_parameters<dim>(ParameterHandler &);

      ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}
