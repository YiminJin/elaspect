#include <elaspect/boundary_composition/interface.h>
#include <elaspect/geometry_model/interface.h>

namespace elaspect
{
  namespace BoundaryComposition
  {

    // ---------------------------- Interface ----------------------------

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


    // ------------------------------ Manager -----------------------------

    template <int dim>
    Manager<dim>::~Manager()
      = default;


    template <int dim>
    void
    Manager<dim>::update()
    {
      for (const auto &p : boundary_composition_objects)
        p->update();
    }


    template <int dim>
    double
    Manager<dim>::boundary_composition(const types::boundary_id boundary_indicator,
                                       const Point<dim> &position,
                                       const unsigned int compositional_field) const
    {
      double composition = 0.0;
      
      for (unsigned int i = 0; i < boundary_composition_objects.size(); ++i)
        composition = model_operators[i](composition,
                                         boundary_composition_objects[i]->boundary_composition(
                                           boundary_indicator,
                                           position,
                                           compositional_field));

      return composition;
    }


    template <int dim>
    const std::set<types::boundary_id> &
    Manager<dim>::get_fixed_composition_boundary_indicators() const
    {
      return fixed_composition_boundary_indicators;
    }


    // -------------------------------- Deal with registering boundary_composition models and automating
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
    Manager<dim>::
    register_boundary_composition (const std::string &name,
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
    Manager<dim>::declare_parameters(ParameterHandler &prm)
    {
      const std::string pattern_of_names
        = std::get<dim>(registered_plugins).get_pattern_of_names();

      prm.enter_subsection("Boundary composition model");
      {
        prm.declare_entry("List of model names", "",
                          Patterns::MultipleSelection(pattern_of_names),
                          "A comma-separated list of boundary composition models that "
                          "will be used to initialize the composition. "
                          "These plugins are loaded in the order given, and modify the "
                          "existing composition field via the operators listed "
                          "in 'List of model operators'.\n\n"
                          "The following boundary composition models are available:\n\n"
                          +
                          std::get<dim>(registered_plugins).get_description_string());

        prm.declare_entry("List of model operators", "add",
                          Patterns::MultipleSelection(Utilities::get_model_operator_options()),
                          "A comma-separated list of operators that "
                          "will be used to append the listed composition models onto "
                          "the previous models. If only one operator is given, "
                          "the same operator is applied to all models.");

        prm.declare_entry ("Fixed composition boundary indicators", "",
                           Patterns::List (Patterns::Anything()),
                           "A comma separated list of names denoting those boundaries "
                           "on which the composition is fixed and described by the "
                           "boundary composition object selected in its own section "
                           "of this input file. All boundary indicators used by the geometry "
"but not explicitly listed here will end up with no-flux "
                           "(insulating) boundary conditions."
                           "\n\n"
                           "The names of the boundaries listed here can either be "
                           "numbers (in which case they correspond to the numerical "
"boundary indicators assigned by the geometry object), or they "
                           "can correspond to any of the symbolic names the geometry object "
                           "may have provided for each part of the boundary. You may want "
                           "to compare this with the documentation of the geometry model you "
                           "use in your model."
                           "\n\n"
                           "This parameter only describes which boundaries have a fixed "
                           "composition, but not what composition should hold on these "
"boundaries. The latter piece of information needs to be "
                           "implemented in a plugin in the BoundaryComposition "
                           "group, unless an existing implementation in this group "
                           "already provides what you want.");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Manager<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary composition model");
      {
        model_names = Utilities::split_string_list(prm.get("List of model names"));

        AssertThrow(Utilities::has_unique_entries(model_names),
                    ExcMessage("The list of strings for the parameter "
                               "'Boundary composition model/List of model names' contains entries more than once. "
                               "This is not allowed. Please check your parameter file."));

        // create operator list
        std::vector<std::string> model_operator_names =
          Utilities::possibly_extend_from_1_to_N(Utilities::split_string_list(prm.get("List of model operators")),
                                                 model_names.size(),
                                                 "List of model operators");
        model_operators = Utilities::create_model_operator_list(model_operator_names);

        try
        {
          const std::vector<types::boundary_id> x_fixed_composition_boundary_indicators
            = this->get_geometry_model().translate_symbolic_boundary_names_to_ids(
                Utilities::split_string_list(prm.get("Fixed composition boundary indicators")));

          fixed_composition_boundary_indicators
            = std::set<types::boundary_id> (x_fixed_composition_boundary_indicators.begin(),
                                            x_fixed_composition_boundary_indicators.end());

          // If model names have been set, but no boundaries on which to use them,
          // ignore the set values, do not create objects that are never used.
          if (fixed_composition_boundary_indicators.size() == 0)
          {
            model_names.clear();
            model_operators.clear();
          }
        }
        catch (const std::string &error)
        {
          AssertThrow (false, ExcMessage ("While parsing the entry <Model settings/Fixed composition "
                                          "boundary indicators>, there was an error. Specifically, "
                                          "the conversion function complained as follows:\n\n"
                                          + error));         
        }
      }
      prm.leave_subsection();

      // go through the list, create objects and let them parse
      // their own parameters
      for (auto &model_name : model_names)
      {
        // create boundary composition objects
        boundary_composition_objects.push_back(std::unique_ptr<Interface<dim>>
                                               (std::get<dim>(registered_plugins)
                                                .create_plugin(model_name,
                                                               "Boundary composition::Model names")));

        if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(boundary_composition_objects.back().get()))
sim->initialize_simulator(this->get_simulator());

        boundary_composition_objects.back()->parse_parameters(prm);
        boundary_composition_objects.back()->initialize();
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
      std::list<internal::Plugins::PluginList<BoundaryComposition::Interface<2>>::PluginInfo> *
      internal::Plugins::PluginList<BoundaryComposition::Interface<2>>::plugins = nullptr;
      template <>
      std::list<internal::Plugins::PluginList<BoundaryComposition::Interface<3>>::PluginInfo> *
      internal::Plugins::PluginList<BoundaryComposition::Interface<3>>::plugins = nullptr;
    }
  }

  namespace BoundaryComposition
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  template class Manager<dim>;

    ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
