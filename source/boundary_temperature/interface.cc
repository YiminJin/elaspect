#include <elaspect/boundary_temperature/interface.h>
#include <elaspect/geometry_model/interface.h>

namespace elaspect
{
  namespace BoundaryTemperature
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
    Interface<dim>::declare_parameters (dealii::ParameterHandler &)
    {}


    template <int dim>
    void
    Interface<dim>::parse_parameters (dealii::ParameterHandler &)
    {}


    // ------------------------------ Manager -----------------------------

    template <int dim>
    Manager<dim>::~Manager()
    {}


    template <int dim>
    void
    Manager<dim>::update ()
    {
      for (unsigned int i=0; i<boundary_temperature_objects.size(); ++i)
        boundary_temperature_objects[i]->update();
    }


    template <int dim>
    double
    Manager<dim>::boundary_temperature (const types::boundary_id boundary_indicator,
                                        const Point<dim> &position) const
    {
      double temperature = 0.0;

      for (unsigned int i=0; i<boundary_temperature_objects.size(); ++i)
        temperature = model_operators[i](temperature,
                                         boundary_temperature_objects[i]->boundary_temperature(boundary_indicator,
                                             position));

      return temperature;
    }


    template <int dim>
    const std::vector<std::string> &
    Manager<dim>::get_active_boundary_temperature_names () const
    {
      return model_names;
    }


    template <int dim>
    const std::vector<std::unique_ptr<Interface<dim> > > &
    Manager<dim>::get_active_boundary_temperature_conditions () const
    {
      return boundary_temperature_objects;
    }


    template <int dim>
    const std::set<types::boundary_id> &
    Manager<dim>::get_fixed_temperature_boundary_indicators() const
    {
      return fixed_temperature_boundary_indicators;
    }


    template <int dim>
    bool
    Manager<dim>::allows_fixed_temperature_on_outflow_boundaries() const
    {
      return allow_fixed_temperature_on_outflow_boundaries;
    }


    // -------------------------------- Deal with registering boundary_temperature models and automating
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
    Manager<dim>::register_boundary_temperature (const std::string &name,
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
      // declare the entry in the parameter file
      prm.enter_subsection ("Boundary temperature model");
      {
        const std::string pattern_of_names
          = std::get<dim>(registered_plugins).get_pattern_of_names ();

        prm.declare_entry("List of model names",
                          "",
                          Patterns::MultipleSelection(pattern_of_names),
                          "A comma-separated list of boundary temperature models that "
                          "will be used to initialize the temperature. "
                          "These plugins are loaded in the order given, and modify the "
                          "existing temperature field via the operators listed "
                          "in 'List of model operators'.\n\n"
                          "The following boundary temperature models are available:\n\n"
                          +
                          std::get<dim>(registered_plugins).get_description_string());

        prm.declare_entry("List of model operators", "add",
                          Patterns::MultipleSelection(Utilities::get_model_operator_options()),
                          "A comma-separated list of operators that "
                          "will be used to append the listed temperature models onto "
                          "the previous models. If only one operator is given, "
                          "the same operator is applied to all models.");

        prm.declare_entry ("Model name", "unspecified",
                           Patterns::Selection (pattern_of_names+"|unspecified"),
                           "Select one of the following models:\n\n"
                           +
                           std::get<dim>(registered_plugins).get_description_string()
                           + "\n\n" +
                           "\\textbf{Warning}: This parameter provides an old and "
                           "deprecated way of specifying "
                           "boundary temperature models and shouldn't be used. "
                           "Please use 'List of model names' instead.");

        prm.declare_entry ("Fixed temperature boundary indicators", "",
                           Patterns::List (Patterns::Anything()),
                           "A comma separated list of names denoting those boundaries "
                           "on which the temperature is fixed and described by the "
                           "boundary temperature object selected in the 'List of model names' "
                           "parameter. All boundary indicators used by the geometry "
                           "but not explicitly listed here will end up with no-flux "
                           "(insulating) boundary conditions, or, if they are listed in the "
                           "'Fixed heat flux boundary indicators', with Neumann boundary "
                           "conditions."
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
                           "temperature, but not what temperature should hold on these "
                           "boundaries. The latter piece of information needs to be "
                           "implemented in a plugin in the BoundaryTemperature "
                           "group, unless an existing implementation in this group "
                           "already provides what you want.");
        prm.declare_entry ("Allow fixed temperature on outflow boundaries", "true",
                           Patterns::Bool (),
                           "When the temperature is fixed on a given boundary as determined "
                           "by the list of 'Fixed temperature boundary indicators', there "
                           "might be parts of the boundary where material flows out and "
                           "one may want to prescribe the temperature only on the parts of "
                           "the boundary where there is inflow. This parameter determines "
                           "if temperatures are only prescribed at these inflow parts of the "
                           "boundary (if false) or everywhere on a given boundary, independent "
                           "of the flow direction (if true)."
                           "Note that in this context, `fixed' refers to the fact that these "
                           "are the boundary indicators where Dirichlet boundary conditions are "
                           "applied, and does not imply that the boundary temperature is "
                           "time-independent. "
                           "\n\n"
                           "Mathematically speaking, the temperature satisfies an "
                           "advection-diffusion equation. For this type of equation, one can "
                           "prescribe the temperature even on outflow boundaries as long as the "
                           "diffusion coefficient is nonzero. This would correspond to the "
                           "``true'' setting of this parameter, which is correspondingly the "
                           "default. In practice, however, this would only make physical sense "
                           "if the diffusion coefficient is actually quite large to prevent "
                           "the creation of a boundary layer. "
                           "In addition, if there is no diffusion, one can only impose "
                           "Dirichlet boundary conditions (i.e., prescribe a fixed temperature "
                           "value at the boundary) at those boundaries where material flows in. "
                           "This would correspond to the ``false'' setting of this parameter.");
      }
      prm.leave_subsection ();

      std::get<dim>(registered_plugins).declare_parameters (prm);
    }


    template <int dim>
    void
    Manager<dim>::parse_parameters (ParameterHandler &prm)
    {
      // find out which plugins are requested and the various other
      // parameters we declare here
      prm.enter_subsection ("Boundary temperature model");
      {
        model_names
          = Utilities::split_string_list(prm.get("List of model names"));

        AssertThrow(Utilities::has_unique_entries(model_names),
                    ExcMessage("The list of strings for the parameter "
                               "'Boundary temperature model/List of model names' contains entries more than once. "
                               "This is not allowed. Please check your parameter file."));

        const std::string model_name = prm.get ("Model name");

        AssertThrow (model_name == "unspecified" || model_names.size() == 0,
                     ExcMessage ("The parameter 'Model name' is only used for reasons"
                                 "of backwards compatibility and can not be used together with "
                                 "the new functionality 'List of model names'. Please add your "
                                 "boundary temperature model to the list instead."));

        if (!(model_name == "unspecified"))
          model_names.push_back(model_name);

        // create operator list
        std::vector<std::string> model_operator_names =
          Utilities::possibly_extend_from_1_to_N (Utilities::split_string_list(prm.get("List of model operators")),
                                                  model_names.size(),
                                                  "List of model operators");
        model_operators = Utilities::create_model_operator_list(model_operator_names);

        try
        {
          const std::vector<types::boundary_id> x_fixed_temperature_boundary_indicators
            = this->get_geometry_model().translate_symbolic_boundary_names_to_ids(Utilities::split_string_list
                                                                                  (prm.get ("Fixed temperature boundary indicators")));
          fixed_temperature_boundary_indicators
            = std::set<types::boundary_id> (x_fixed_temperature_boundary_indicators.begin(),
                                            x_fixed_temperature_boundary_indicators.end());

          // If model names have been set, but no boundaries on which to use them,
          // ignore the set values, do not create objects that are never used.
          if (fixed_temperature_boundary_indicators.size() == 0)
          {
            model_names.clear();
            model_operators.clear();
          }
      }
        catch (const std::string &error)
        {
          AssertThrow (false, ExcMessage ("While parsing the entry <Model settings/Fixed temperature "
                                          "boundary indicators>, there was an error. Specifically, "
                                          "the conversion function complained as follows: "
                                          + error));
        }

        allow_fixed_temperature_on_outflow_boundaries = prm.get_bool ("Allow fixed temperature on outflow boundaries");
      }
      prm.leave_subsection ();

      // go through the list, create objects and let them parse
      // their own parameters
      for (auto &model_name : model_names)
      {
        // create boundary temperature objects
        boundary_temperature_objects.push_back (std::unique_ptr<Interface<dim> >
                                                (std::get<dim>(registered_plugins)
                                                 .create_plugin (model_name,
                                                                 "Boundary temperature::Model names")));

        if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(boundary_temperature_objects.back().get()))
          sim->initialize_simulator (this->get_simulator());

        boundary_temperature_objects.back()->parse_parameters (prm);
        boundary_temperature_objects.back()->initialize ();
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
      std::list<internal::Plugins::PluginList<BoundaryTemperature::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<BoundaryTemperature::Interface<2> >::plugins = nullptr;
      template <>
      std::list<internal::Plugins::PluginList<BoundaryTemperature::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<BoundaryTemperature::Interface<3> >::plugins = nullptr;
    }
  }

  namespace BoundaryTemperature
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  template class Manager<dim>;

    ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
