#include <elaspect/geometry_model/sphere.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

namespace elaspect
{
  namespace GeometryModel
  {
    template <int dim>
    Sphere<dim>::Sphere ()
      : spherical_manifold ()
    {}


    template <int dim>
    void
    Sphere<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
      GridGenerator::hyper_ball_balanced (coarse_grid,
                                          Point<dim>(),
                                          R);

      coarse_grid.set_manifold(0,spherical_manifold);
      coarse_grid.set_all_manifold_ids_on_boundary(0);
    }


    template <int dim>
    std::set<types::boundary_id>
    Sphere<dim>::
    get_used_boundary_indicators () const
    {
      const types::boundary_id s[] = { 0 };
      return std::set<types::boundary_id>(&s[0],
                                          &s[sizeof(s)/sizeof(s[0])]);
    }


    template <int dim>
    std::map<std::string,types::boundary_id>
    Sphere<dim>::
    get_symbolic_boundary_names_map () const
    {
      static const std::pair<std::string,types::boundary_id> mapping("top", 0);
      return std::map<std::string,types::boundary_id> (&mapping,
                                                       &mapping+1);
    }


    template <int dim>
    double
    Sphere<dim>::depth(const Point<dim> &position) const
    {
      return std::min (std::max (R-position.norm(), 0.), maximal_depth());
    }


    template <int dim>
    double
    Sphere<dim>::maximal_depth() const
    {
      return R;
    }


    template <int dim>
    double
    Sphere<dim>::height_above_reference_surface(const Point<dim> &position) const
    {
      return position.norm() - radius();
    }


    template <int dim>
    elaspect::Utilities::Coordinates::CoordinateSystem
    Sphere<dim>::natural_coordinate_system() const
    {
      return elaspect::Utilities::Coordinates::CoordinateSystem::spherical;
    }


    template <int dim>
    std::array<double,dim>
    Sphere<dim>::cartesian_to_natural_coordinates(const Point<dim> &position) const
    {
      return Utilities::Coordinates::cartesian_to_spherical_coordinates<dim>(position);
    }


    template <int dim>
    Point<dim>
    Sphere<dim>::natural_to_cartesian_coordinates(const std::array<double,dim> &position) const
    {
      return Utilities::Coordinates::spherical_to_cartesian_coordinates<dim>(position);
    }


    template <int dim>
    double Sphere<dim>::radius () const
    {
      return R;
    }


    template <int dim>
    void
    Sphere<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Sphere");
        {
          prm.declare_entry ("Radius", "6371000.",
                             Patterns::Double (0.),
                             "Radius of the sphere. Units: \\si{\\meter}.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Sphere<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Sphere");
        {
          R = prm.get_double ("Radius");
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
  namespace GeometryModel
  {
    ELASPECT_REGISTER_GEOMETRY_MODEL(Sphere,
                                   "sphere",
                                   "")
  }
}
