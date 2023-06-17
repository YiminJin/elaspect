#include <elaspect/geometry_model/box.h>
#include <elaspect/initial_topography/zero_topography.h>
#include <elaspect/simulator_signals.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

namespace elaspect
{
  namespace GeometryModel
  {
    template <int dim>
    void Box<dim>::initialize()
    {
      const InitialTopography::Interface<dim> &topo_model = this->get_initial_topography();

      if (!Plugins::plugin_type_matches<InitialTopography::ZeroTopography<dim>>(topo_model))
      {
        this->get_signals().pre_set_initial_state.connect(
          [&](parallel::distributed::Triangulation<dim> &tria)
        {
          this->add_topography(topo_model, tria);
        });
      }
    }


    template <int dim>
    void
    Box<dim>::
    add_topography(const InitialTopography::Interface<dim> &topo_model,
                   parallel::distributed::Triangulation<dim> &grid) const
    {
      auto add_topo = [&] (const Point<dim> &p) -> Point<dim>
      {
        Point<dim-1> surface_point;
        for (unsigned int d = 0; d < dim-1; ++d)
          surface_point[d] = p[d];

        const double topo = topo_model.value(surface_point);

        const double ztopo = (p[dim-1] - box_origin[dim-1]) / extents[dim-1] * topo;

        Point<dim> result = p;
        result[dim-1] += ztopo;

        return result;
      };

      GridTools::transform(add_topo, grid);

      this->get_pcout() << "   Added initial topography to grid." << std::endl << std::endl;
    }


    template <int dim>
    void
    Box<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
      std::vector<unsigned int> rep_vec(repetitions, repetitions+dim);
      GridGenerator::subdivided_hyper_rectangle (coarse_grid,
                                                 rep_vec,
                                                 box_origin,
                                                 box_origin+extents,
                                                 true);
    }


    template <int dim>
    std::set<types::boundary_id>
    Box<dim>::
    get_used_boundary_indicators () const
    {
      // boundary indicators are zero through 2*dim-1
      std::set<types::boundary_id> s;
      for (unsigned int i=0; i<2*dim; ++i)
        s.insert (i);
      return s;
    }


    template <int dim>
    std::map<std::string,types::boundary_id>
    Box<dim>::
    get_symbolic_boundary_names_map () const
    {
      switch (dim)
        {
          case 2:
          {
            static const std::pair<std::string,types::boundary_id> mapping[]
              = { std::pair<std::string,types::boundary_id>("left",   0),
                  std::pair<std::string,types::boundary_id>("right",  1),
                  std::pair<std::string,types::boundary_id>("bottom", 2),
                  std::pair<std::string,types::boundary_id>("top",    3)
                };

            return std::map<std::string,types::boundary_id> (&mapping[0],
                                                             &mapping[sizeof(mapping)/sizeof(mapping[0])]);
          }

          case 3:
          {
            static const std::pair<std::string,types::boundary_id> mapping[]
              = { std::pair<std::string,types::boundary_id>("left",   0),
                  std::pair<std::string,types::boundary_id>("right",  1),
                  std::pair<std::string,types::boundary_id>("front",  2),
                  std::pair<std::string,types::boundary_id>("back",   3),
                  std::pair<std::string,types::boundary_id>("bottom", 4),
                  std::pair<std::string,types::boundary_id>("top",    5)
                };

            return std::map<std::string,types::boundary_id> (&mapping[0],
                                                             &mapping[sizeof(mapping)/sizeof(mapping[0])]);
          }
        }

      Assert (false, ExcNotImplemented());
      return std::map<std::string,types::boundary_id>();
    }


    template <int dim>
    double
    Box<dim>::depth(const Point<dim> &position) const
    { 
      const double d = extents[dim-1] - (position(dim-1)-box_origin[dim-1]);
      return std::min (std::max (d, 0.), maximal_depth());
    }


    template <int dim>
    double
    Box<dim>::maximal_depth() const
    {
      return extents[dim-1];
    }


    template <int dim>
    double
    Box<dim>::height_above_reference_surface(const Point<dim> &position) const
    {
      return (position[dim-1] - box_origin[dim-1]) - extents[dim-1];
    }


    template <int dim>
    std::array<double,dim>
    Box<dim>::cartesian_to_natural_coordinates(const Point<dim> &position_point) const
    {
      std::array<double,dim> position_array;
      for (unsigned int i = 0; i < dim; i++)
        position_array[i] = position_point(i);

      return position_array;
    }


    template <int dim>
    elaspect::Utilities::Coordinates::CoordinateSystem
    Box<dim>::natural_coordinate_system() const
    {
      return elaspect::Utilities::Coordinates::CoordinateSystem::cartesian;
    }


    template <int dim>
    Point<dim>
    Box<dim>::natural_to_cartesian_coordinates(const std::array<double,dim> &position_tensor) const
    {
      Point<dim> position_point;
      for (unsigned int i = 0; i < dim; i++)
        position_point[i] = position_tensor[i];

      return position_point;
    }


    template <int dim>
    void
    Box<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Box");
        {
          prm.declare_entry ("X extent", "1.",
                             Patterns::Double (0.),
                             "Extent of the box in x-direction. Units: \\si{\\meter}.");
          prm.declare_entry ("Y extent", "1.",
                             Patterns::Double (0.),
                             "Extent of the box in y-direction. Units: \\si{\\meter}.");
          prm.declare_entry ("Z extent", "1.",
                             Patterns::Double (0.),
                             "Extent of the box in z-direction. This value is ignored "
                             "if the simulation is in 2d. Units: \\si{\\meter}.");

          prm.declare_entry ("Box origin X coordinate", "0.",
                             Patterns::Double (),
                             "X coordinate of box origin. Units: \\si{\\meter}.");
          prm.declare_entry ("Box origin Y coordinate", "0.",
                             Patterns::Double (),
                             "Y coordinate of box origin. Units: \\si{\\meter}.");
          prm.declare_entry ("Box origin Z coordinate", "0.",
                             Patterns::Double (),
                             "Z coordinate of box origin. This value is ignored "
                             "if the simulation is in 2d. Units: \\si{\\meter}.");

          prm.declare_entry ("X repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in X direction.");
          prm.declare_entry ("Y repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in Y direction.");
          prm.declare_entry ("Z repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in Z direction.");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Box<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Box");
        {
          box_origin[0] = prm.get_double ("Box origin X coordinate");
          extents[0] = prm.get_double ("X extent");
          repetitions[0] = prm.get_integer ("X repetitions");

          if (dim >= 2)
            {
              box_origin[1] = prm.get_double ("Box origin Y coordinate");
              extents[1] = prm.get_double ("Y extent");
              repetitions[1] = prm.get_integer ("Y repetitions");
            }

          if (dim >= 3)
            {
              // Use dim-1 instead of 2 to avoid compiler warning in 2d:
              box_origin[dim-1] = prm.get_double ("Box origin Z coordinate");
              extents[dim-1] = prm.get_double ("Z extent");
              repetitions[dim-1] = prm.get_integer ("Z repetitions");
            }
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
    ELASPECT_REGISTER_GEOMETRY_MODEL(Box,
                                 "box",
                                 "")
  }
}
