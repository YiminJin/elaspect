#include <elaspect/initial_composition/function.h>
#include <elaspect/geometry_model/interface.h>

namespace elaspect
{
  namespace InitialComposition
  {
    template <int dim>
    double
    Function<dim>::
    initial_composition (const Point<dim> &position, 
                         const unsigned int n_comp) const
    {
      Utilities::NaturalCoordinate<dim> point =
        this->get_geometry_model().cartesian_to_other_coordinates(position, coordinate_system);

      return function->value(Utilities::convert_array_to_point<dim>(point.get_coordinates()), n_comp);
    }


    template <int dim>
    void
    Function<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Function");
        {
          /**
           * Choose the coordinates to evaluate the maximum refinement level
           * function. The function can be declared in dependence of depth,
           * cartesian coordinates or spherical coordinates. Note that the order
           * of spherical coordinates is r,phi,theta and not r,theta,phi, since
           * this allows for dimension independent expressions.
           */
          prm.declare_entry ("Coordinate system", "cartesian",
                             Patterns::Selection ("cartesian|spherical|depth"),
                             "A selection that determines the assumed coordinate "
                             "system for the function variables. Allowed values "
                             "are `cartesian', `spherical', and `depth'. `spherical' coordinates "
                             "are interpreted as r,phi or r,phi,theta in 2D/3D "
                             "respectively with theta being the polar angle. `depth' "
                             "will create a function, in which only the first "
                             "parameter is non-zero, which is interpreted to "
                             "be the depth of the point.");

          Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Function<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Function");
        {
          coordinate_system = Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));
        }

        try
        {
          function = std::make_unique<Functions::ParsedFunction<dim> >(this->n_compositional_fields());
          function->parse_parameters (prm);
        }
        catch (...)
        {
          std::cerr << "ERROR: FunctionParser failed to parse\n"
                    << "\t'Initial composition model.Function'\n"
                    << "with expression\n"
                    << "\t'" << prm.get("Function expression") << "'\n"
                    << "More information about the cause of the parse error \n"
                    << "is shown below.\n";
          throw;
        }

        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace InitialComposition
  {
    ELASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(Function,
                                            "function",
                                            "")
  }
}
