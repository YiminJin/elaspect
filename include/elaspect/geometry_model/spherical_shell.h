#ifndef _elaspect_geometry_model_spherical_shell_h
#define _elaspect_geometry_model_spherical_shell_h

#include <elaspect/geometry_model/interface.h>
#include <elaspect/simulator_access.h>

#include <deal.II/grid/manifold_lib.h>

namespace elaspect
{
  namespace GeometryModel
  {
    using namespace dealii;

    template <int dim>
    class SphericalShell : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        SphericalShell ();

        void create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const override;

        std::set<types::boundary_id>
        get_used_boundary_indicators () const override;

        std::map<std::string,types::boundary_id>
        get_symbolic_boundary_names_map () const override;

        double depth(const Point<dim> &position) const override;

        double maximal_depth() const override;

        double height_above_reference_surface(const Point<dim> &position) const override;

        elaspect::Utilities::Coordinates::CoordinateSystem 
        natural_coordinate_system() const override;

        std::array<double,dim> 
        cartesian_to_natural_coordinates(const Point<dim> &position) const override;

        Point<dim> 
        natural_to_cartesian_coordinates(const std::array<double,dim> &position) const override;

        double inner_radius () const;

        double outer_radius () const;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        enum CustomMeshRadialSubdivision
        {
          none,
          list,
          slices
        } custom_mesh;

        int initial_lateral_refinement;

        unsigned int n_slices;

        std::vector<double> R_values_list;

        double R0, R1;

        double phi;

        int n_cells_along_circumference;

        const SphericalManifold<dim> spherical_manifold;

        void set_manifold_ids (parallel::distributed::Triangulation<dim> &triangulation) const;
    };
  }
}

#endif
