#ifndef _elaspect_geometry_model_sphere_h
#define _elaspect_geometry_model_sphere_h

#include <elaspect/geometry_model/interface.h>
#include <elaspect/simulator_access.h>

#include <deal.II/grid/manifold_lib.h>

namespace elaspect
{
  namespace GeometryModel
  {
    using namespace dealii;

    template <int dim>
    class Sphere : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        Sphere ();

        void create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const override;

        std::set<types::boundary_id>
        get_used_boundary_indicators () const override;

        std::map<std::string,types::boundary_id>
        get_symbolic_boundary_names_map () const override;

        double depth(const Point<dim> &position) const override;

        double maximal_depth() const override;

        double height_above_reference_surface(const Point<dim> &position) const override;

        elaspect::Utilities::Coordinates::CoordinateSystem
        natural_coordinate_system () const override;

        std::array<double,dim> 
        cartesian_to_natural_coordinates (const Point<dim> &position) const override;

        Point<dim>
        natural_to_cartesian_coordinates (const std::array<double,dim> &position) const override;

        double radius () const;

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        double R;

        const SphericalManifold<dim> spherical_manifold;
    };
  }
}

#endif
