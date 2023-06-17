#include <elaspect/boundary_heat_flux/function.h>
#include <elaspect/geometry_model/interface.h>

namespace elaspect
{
  namespace BoundaryHeatFlux
  {
    template <int dim>
    Tensor<1,dim>
    Function<dim>::heat_flux(const types::boundary_id /*boundary_id*/,
                             const Point<dim>    &position,
                             const Tensor<1,dim> &normal_vector) const
    {
      Tensor<1,dim> flux = normal_vector;
      if (coordinate_system == Utilities::Coordinates::cartesian)
      {
        flux *= boundary_heat_flux_function.value(position);
      }
      else if (coordinate_system == Utilities::Coordinates::spherical)
      {
        const std::array<double,dim> spherical_coordinates =
          Utilities::Coordinates::cartesian_to_spherical_coordinates(position);

        Point<dim> point;
        for (unsigned int d = 0; d < dim; ++d)
          point[d] = spherical_coordinates[d];

        flux *= boundary_heat_flux_function.value(point);
      }
      else if (coordinate_system == Utilities::Coordinates::depth)
      {
        const double depth = this->get_geometry_model().depth(position);
        Point<dim> point;
        point(0) = depth;

        flux *= boundary_heat_flux_function.value(point);
      }
      else
      {
        Assert(false, ExcNotImplemented());
      }

      return flux;
    }


    template <int dim>
    void
    Function<dim>::update()
    {
      // we get time passed as seconds (always) but may want
      // to reinterpret it in years
      if (this->convert_output_to_years())
        boundary_heat_flux_function.set_time (this->get_time() / year_in_seconds);
      else
        boundary_heat_flux_function.set_time (this->get_time());
    }


    template <int dim>
    void
    Function<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary heat flux model");
      {
        prm.enter_subsection("Function");
        {
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
      prm.enter_subsection("Boundary heat flux model");
      {
        prm.enter_subsection("Function");
        {
          coordinate_system = Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));
        }
        try
          {
            boundary_heat_flux_function.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Boundary heat flux model.Function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'"
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
  namespace BoundaryHeatFlux
  {
    ELASPECT_REGISTER_BOUNDARY_HEAT_FLUX_MODEL(Function,
                                           "function",
                                           "Implementation of a model in which the boundary "
                                           "heat flux is given in terms of an explicit formula "
                                           "that is elaborated in the parameters in section "
                                           "``Boundary heat flux model|Function''. The format of these "
                                           "functions follows the syntax understood by the "
                                           "muparser library, see Section~\\ref{sec:muparser-format}."
                                           "\n\n"
                                           "The formula you describe in the mentioned "
                                           "section is a scalar value for the heat flux that is assumed "
                                           "to be the flux normal to the boundary, and that has the unit "
                                           "W/(m$^2$) (in 3d) or W/m (in 2d). Negative fluxes are "
                                           "interpreted as the flow of heat into the domain, and positive "
                                           "fluxes are interpreted as heat flowing out of the domain."
                                           "\n\n"
                                           "The symbol $t$ indicating time that "
                                           "may appear in the formulas for the prescribed "
                                           "heat flux is interpreted as having units "
                                           "seconds unless the global parameter ``Use "
                                           "years in output instead of seconds'' has "
                                           "been set.")
  }
}
