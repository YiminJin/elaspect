#include <elaspect/geometry_model/interface.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

namespace elaspect
{
  namespace GeometryModel
  {
    template <int dim>
    class SlopeStability : public Interface<dim>
    {
      public:
        void 
        create_coarse_mesh(parallel::distributed::Triangulation<dim> &coarse_grid) const override;

        double depth(const Point<dim> &position) const override;

        double maximal_depth() const override;

        double height_above_reference_surface(const Point<dim> &position) const override;

        std::set<types::boundary_id>
        get_used_boundary_indicators() const override;

        std::map<std::string, types::boundary_id>
        get_symbolic_boundary_names_map() const override;

        elaspect::Utilities::Coordinates::CoordinateSystem natural_coordinate_system() const override;

        std::array<double,dim> 
        cartesian_to_natural_coordinates(const Point<dim> &position) const override;

        Point<dim> 
        natural_to_cartesian_coordinates(const std::array<double,dim> &position) const override;
    };


    template <int dim>
    void
    SlopeStability<dim>::
    create_coarse_mesh(parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
      unsigned int repetitions[dim];
      
      // create the slope
      Triangulation<dim> slope;
      Point<dim> p1, p2;
      p2[0]              = 20;
      p2[dim-1]          = 10;
      repetitions[0]     = 4;
      repetitions[dim-1] = 2;
      if (dim == 3)
      {
        p2[1]            = 50;
        repetitions[1]   = 5;
      }

      GridGenerator::subdivided_hyper_rectangle(slope, 
                                                std::vector<unsigned int>(repetitions,
                                                                          repetitions+dim), 
                                                p1, 
                                                p2, 
                                                false);

      auto trans = [&] (const Point<dim> &p) -> Point<dim>
      {
        Point<dim> p_new(p);
        p_new[0] = p[0] * (-0.1 * p[dim-1] + 2);
        return p_new;
      };
      GridTools::transform(trans, slope);

      // create the base
      Triangulation<dim> base;
      Point<dim> p3, p4;
      p3[dim-1]          = -5;
      p4[0]              = 60;
      repetitions[0]     = 6;
      repetitions[dim-1] = 1;
      if (dim == 3)
      {
        p4[1]            = 50;
        repetitions[1]   = 5;
      }

      GridGenerator::subdivided_hyper_rectangle(base, 
                                                std::vector<unsigned int>(repetitions,
                                                                          repetitions+dim), 
                                                p3, 
                                                p4, 
                                                false);

      // merge the two triangulations
      const double tol = 1e-4;
      GridGenerator::merge_triangulations(base, slope, coarse_grid, tol);

      for (const auto &cell : coarse_grid.active_cell_iterators())
      {
        if (std::fabs(cell->face(0)->center()(0)) < tol)
          cell->face(0)->set_boundary_id(10);
        if (std::fabs(cell->face(1)->center()(0) - 60) < tol)
          cell->face(1)->set_boundary_id(11);

        if (dim == 2)
        {
          if (std::fabs(cell->face(2)->center()(1) + 5) < 1e-4)
            cell->face(2)->set_boundary_id(12);
          if (std::fabs(cell->face(3)->center()(1) - 10) < 1e-4)
            cell->face(3)->set_boundary_id(13);
        }
        else
        {
          if (std::fabs(cell->face(2)->center()(1)) < tol)
            cell->face(2)->set_boundary_id(12);
          if (std::fabs(cell->face(3)->center()(1) - 50) < tol)
            cell->face(3)->set_boundary_id(13);
          if (std::fabs(cell->face(4)->center()(2) + 5) < tol)
            cell->face(4)->set_boundary_id(14);
          if (std::fabs(cell->face(5)->center()(2) - 10) < tol)
            cell->face(5)->set_boundary_id(15);
        }
      }
    }


    template <int dim>
    std::set<types::boundary_id>
    SlopeStability<dim>::get_used_boundary_indicators() const
    {
      std::set<types::boundary_id> s;
      // default boundary indicator is 0
      s.insert(0);
      for (unsigned int i = 0; i < 2*dim; ++i)
        s.insert(10 + i);
      return s;
    }


    template <int dim>
    std::map<std::string, types::boundary_id>
    SlopeStability<dim>::get_symbolic_boundary_names_map() const
    {
      switch (dim)
      {
        case 2:
        {
          static const std::pair<std::string, types::boundary_id> mapping[]
            = { std::pair<std::string, types::boundary_id>("left",   10),
                std::pair<std::string, types::boundary_id>("right",  11),
                std::pair<std::string, types::boundary_id>("bottom", 12),
                std::pair<std::string, types::boundary_id>("top",    13)
              };

          return std::map<std::string, types::boundary_id>(
            &mapping[0], &mapping[sizeof(mapping)/sizeof(mapping[0])]);
        }
        case 3:
        {
          static const std::pair<std::string, types::boundary_id> mapping[]
            = { std::pair<std::string, types::boundary_id>("left",   10),
                std::pair<std::string, types::boundary_id>("right",  11),
                std::pair<std::string, types::boundary_id>("front",  12),
                std::pair<std::string, types::boundary_id>("back",   13),
                std::pair<std::string, types::boundary_id>("bottom", 14),
                std::pair<std::string, types::boundary_id>("top",    15)
              };

          return std::map<std::string, types::boundary_id>(
            &mapping[0], &mapping[sizeof(mapping)/sizeof(mapping[0])]);
        }
      }

      Assert(false, ExcNotImplemented());
      return std::map<std::string, types::boundary_id>();
    }


    template <int dim>
    double SlopeStability<dim>::depth(const Point<dim> &position) const
    {
      return std::min(std::max(-position[dim-1], 0.), maximal_depth());
    }


    template <int dim>
    double SlopeStability<dim>::maximal_depth() const
    {
      return 5.;
    }


    template <int dim>
    double
    SlopeStability<dim>::height_above_reference_surface(const Point<dim> &position) const
    {
      return 10 - position[dim-1];
    }


    template <int dim>
    std::array<double, dim>
    SlopeStability<dim>::
    cartesian_to_natural_coordinates(const Point<dim> &position_point) const
    {
      std::array<double, dim> position_array;
      for (unsigned int i = 0; i < dim; i++)
        position_array[i] = position_point(i);

      return position_array;
    }


    template <int dim>
    elaspect::Utilities::Coordinates::CoordinateSystem
    SlopeStability<dim>::natural_coordinate_system() const
    {
      return elaspect::Utilities::Coordinates::CoordinateSystem::cartesian;
    }


    template <int dim>
    Point<dim>
    SlopeStability<dim>::
    natural_to_cartesian_coordinates(const std::array<double, dim> &position_tensor) const
    {
      Point<dim> position_point;
      for (unsigned int i = 0; i < dim; i++)
        position_point[i] = position_tensor[i];

      return position_point;
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace GeometryModel
  {
    ELASPECT_REGISTER_GEOMETRY_MODEL(SlopeStability,
                                 "slope stability",
                                 "")
  }
}
