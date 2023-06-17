#include <elaspect/geometry_model/interface.h>
#include <elaspect/postprocess/interface.h>
#include <elaspect/simulator_access.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/grid/grid_generator.h>

namespace elaspect
{
  namespace GeometryModel
  {
    template <int dim>
    class StripFooting : public Interface<dim>
    {
      public:
        void create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const override;

        double depth (const Point<dim> &position) const override;

        double maximal_depth () const override;

        std::set<types::boundary_id>
        get_used_boundary_indicators () const override;

        std::map<std::string,types::boundary_id>
        get_symbolic_boundary_names_map () const override;

        elaspect::Utilities::Coordinates::CoordinateSystem natural_coordinate_system() const override;

        std::array<double,dim> cartesian_to_natural_coordinates(const Point<dim> &position) const override;

        Point<dim> natural_to_cartesian_coordinates(const std::array<double,dim> &position) const override;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        Point<dim> extents;

        unsigned int repetitions[dim];

        double strip_half_width;
    };


    template <int dim>
    void
    StripFooting<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_mesh) const
    {
      std::vector<unsigned int> rep_vec(repetitions, repetitions+dim);
      const Point<dim> box_origin;
      GridGenerator::subdivided_hyper_rectangle (coarse_mesh, 
                                                 rep_vec, 
                                                 box_origin, 
                                                 box_origin + extents, 
                                                 true);

      // colorize the strip footing area
      const types::boundary_id top_boundary_id = (dim == 2 ? 3 : 5);
      for (const auto &cell : coarse_mesh.active_cell_iterators())
        if (cell->at_boundary(top_boundary_id))
          if (cell->center()(0) <= strip_half_width)
            cell->face(top_boundary_id)->set_boundary_id(99);
    }


    template <int dim>
    std::set<types::boundary_id>
    StripFooting<dim>::get_used_boundary_indicators () const
    {
      std::set<types::boundary_id> s;
      for (unsigned int i = 0; i < 2 * dim; ++i)
        s.insert(i);
      s.insert(99);
      return s;
    }


    template <int dim>
    std::map<std::string, types::boundary_id>
    StripFooting<dim>::get_symbolic_boundary_names_map () const
    {
      switch (dim)
      {
        case 2:
        {
          static const std::pair<std::string,types::boundary_id> mapping[]
            = { std::pair<std::string,types::boundary_id>("left",   0),
                std::pair<std::string,types::boundary_id>("right",  1),
                std::pair<std::string,types::boundary_id>("bottom", 2),
                std::pair<std::string,types::boundary_id>("top",    3),
                std::pair<std::string,types::boundary_id>("strip", 99)
              };

          return std::map<std::string,types::boundary_id> (&mapping[0],
                                                           &mapping[5]);
        }
        case 3:
        {
          static const std::pair<std::string,types::boundary_id> mapping[]
            = { std::pair<std::string,types::boundary_id>("left",   0),
                std::pair<std::string,types::boundary_id>("right",  1),
                std::pair<std::string,types::boundary_id>("front",  2),
                std::pair<std::string,types::boundary_id>("back",   3),
                std::pair<std::string,types::boundary_id>("bottom", 4),
                std::pair<std::string,types::boundary_id>("top",    5),
                std::pair<std::string,types::boundary_id>("strip", 99)
              };

          return std::map<std::string,types::boundary_id> (&mapping[0],
                                                           &mapping[7]);
        }
      }

      Assert (false, ExcNotImplemented());
      return std::map<std::string,types::boundary_id>();
    }


    template <int dim>
    double
    StripFooting<dim>::depth (const Point<dim> &position) const
    {
      const double d = extents[dim-1] - (position(dim-1));
      return std::min (std::max (d, 0.), maximal_depth());
    }


    template <int dim>
    double
    StripFooting<dim>::maximal_depth () const
    {
      return extents[dim-1];
    }


    template <int dim>
    std::array<double,dim>
    StripFooting<dim>::cartesian_to_natural_coordinates(const Point<dim> &position_point) const
    {
      std::array<double,dim> position_array;
      for (unsigned int i = 0; i < dim; i++)
        position_array[i] = position_point(i);

      return position_array;
    }


    template <int dim>
    elaspect::Utilities::Coordinates::CoordinateSystem
    StripFooting<dim>::natural_coordinate_system() const
    {
      return elaspect::Utilities::Coordinates::CoordinateSystem::cartesian;
    }


    template <int dim>
    Point<dim>
    StripFooting<dim>::natural_to_cartesian_coordinates(const std::array<double,dim> &position_tensor) const
    {
      Point<dim> position_point;
      for (unsigned int i = 0; i < dim; i++)
        position_point[i] = position_tensor[i];

      return position_point;
    }


    template <int dim>
    void
    StripFooting<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Strip footing");
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

          prm.declare_entry ("X repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in X direction.");
          prm.declare_entry ("Y repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in Y direction.");
          prm.declare_entry ("Z repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in Z direction.");

          prm.declare_entry ("Strip half width", "1",
                             Patterns::Double (0),
                             "");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    StripFooting<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Strip footing");
        {
          extents[0] = prm.get_double ("X extent");
          extents[1] = prm.get_double ("Y extent");

          repetitions[0] = prm.get_integer ("X repetitions");
          repetitions[1] = prm.get_integer ("Y repetitions");

          if (dim == 3)
          {
            extents[dim-1] = prm.get_double ("Z extent");
            repetitions[dim-1] = prm.get_integer ("Z repetitions");
          }

          strip_half_width = prm.get_double ("Strip half width");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }


  namespace Postprocess
  {
    template <int dim>
    class BearingCapacity : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        std::pair<std::string,std::string>
        execute(TableHandler &statistics) override;
    };


    template <int dim>
    std::pair<std::string,std::string>
    BearingCapacity<dim>::execute(TableHandler &statistics)
    {
      const QGauss<dim-1> face_quadrature_formula(this->get_parameters().n_gaussian_points);
      FEFaceValues<dim> fe_face_values(this->get_mapping(), this->get_fe(),
                                       face_quadrature_formula,
                                       update_values | update_JxW_values);

      const unsigned int n_q_points = face_quadrature_formula.size();
      std::vector<SymmetricTensor<2,dim>> stress_values(n_q_points);

      // hardcoding here
      const unsigned int top_face_no = (dim == 2 ? 3 : 5);

      double force_local = 0;
      double area_local = 0;
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          if (cell->at_boundary(top_face_no) &&
              cell->face(top_face_no)->center()(0) <= 1)
          {
            fe_face_values.reinit(cell, top_face_no);
            fe_face_values[this->introspection().extractors.stress].get_function_values(
              this->get_solution(), stress_values);

            for (unsigned int q = 0; q < n_q_points; ++q)
            {
              force_local += stress_values[q][dim-1][dim-1] * fe_face_values.JxW(q);
              area_local += fe_face_values.JxW(q);
            }
          }

      const double force_global = Utilities::MPI::sum(force_local, this->get_mpi_communicator());
      const double area_global  = Utilities::MPI::sum(area_local, this->get_mpi_communicator());

      const double b_c = force_global / area_global;

      const std::string column_name("Bearing capacity");
      statistics.add_value(column_name, b_c);
      statistics.set_precision(column_name, 8);
      statistics.set_scientific(column_name, true);

      std::ostringstream output;
      output.precision(4);
      output << b_c;

      return std::pair<std::string, std::string>("Bearing capacity:", output.str());
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace GeometryModel
  {
    ELASPECT_REGISTER_GEOMETRY_MODEL(StripFooting,
                                 "strip footing",
                                 "")
  }

  namespace Postprocess
  {
    ELASPECT_REGISTER_POSTPROCESSOR(BearingCapacity,
                                "bearing capacity",
                                "")
  }
}
