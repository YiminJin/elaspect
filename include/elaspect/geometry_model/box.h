/*
  Copyright (C) 2023 by Yimin Jin.

  This file is part of elASPECT.

  elASPECT is modified from the free software ASPECT; you can 
  redistribute it and/or modify it under the terms of the GNU 
  General Public License as published by the Free Software 
  Foundation; either version 2, or (at your option) any later 
  version.

  elASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with elASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef _elaspect_geometry_model_box_h
#define _elaspect_geometry_model_box_h

#include <elaspect/geometry_model/interface.h>
#include <elaspect/initial_topography/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace GeometryModel
  {
    using namespace dealii;

    template <int dim>
    class Box : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        void initialize() override;

        void create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const override;

        double depth (const Point<dim> &position) const override;

        double maximal_depth () const override;

        double height_above_reference_surface(const Point<dim> &position) const override;

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
        void add_topography (const InitialTopography::Interface<dim> &topo_model,
                             parallel::distributed::Triangulation<dim> &grid) const;

        Point<dim> extents;

        Point<dim> box_origin;

        unsigned int repetitions[dim];
    };
  }
}

#endif
