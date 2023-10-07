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


#ifndef _elaspect_geometry_model_interface_h
#define _elaspect_geometry_model_interface_h

#include <elaspect/plugins.h>
#include <elaspect/utilities.h>

#include <deal.II/distributed/tria.h>

namespace elaspect
{
  namespace GeometryModel
  {
    using namespace dealii;

    template <int dim>
    class Interface
    {
      public:
        virtual ~Interface ();

        virtual void initialize ();

        virtual
        void create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const = 0;
        
        virtual
        std::set<types::boundary_id>
        get_used_boundary_indicators () const = 0;

        virtual
        std::map<std::string,types::boundary_id>
        get_symbolic_boundary_names_map () const;

        types::boundary_id
        translate_symbolic_boundary_name_to_id (const std::string &name) const;
        
        std::vector<types::boundary_id>
        translate_symbolic_boundary_names_to_ids (const std::vector<std::string> &names) const;

        std::string
        translate_id_to_symbol_name (const types::boundary_id boundary_id) const;

        virtual
        double depth(const Point<dim> &position) const = 0;

        virtual
        double maximal_depth() const = 0;

        virtual
        double height_above_reference_surface(const Point<dim> &position) const = 0;

        Utilities::NaturalCoordinate<dim>
        cartesian_to_other_coordinates(const Point<dim> &position,
                                       const Utilities::Coordinates::CoordinateSystem &coordinate_system) const;

        virtual
        elaspect::Utilities::Coordinates::CoordinateSystem 
        natural_coordinate_system() const = 0;

        virtual
        std::array<double,dim> 
        cartesian_to_natural_coordinates(const Point<dim> &position) const;

        virtual
        Point<dim> 
        natural_to_cartesian_coordinates(const std::array<double,dim> &position) const;

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);
    };

    template <int dim>
    void
    register_geometry_model (const std::string &name,
                             const std::string &description,
                             void (*declare_parameters_function) (ParameterHandler &),
                             Interface<dim> *(*factory_function) ());

    template <int dim>
    Interface<dim> *
    create_geometry_model (ParameterHandler &prm);

    template <int dim>
    void
    declare_parameters (ParameterHandler &prm);


#define ELASPECT_REGISTER_GEOMETRY_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ELASPECT_REGISTER_GEOMETRY_MODEL_ ## classname \
  { \
    elaspect::internal::Plugins::RegisterHelper<elaspect::GeometryModel::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&elaspect::GeometryModel::register_geometry_model<2>, \
                                name, description); \
    elaspect::internal::Plugins::RegisterHelper<elaspect::GeometryModel::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&elaspect::GeometryModel::register_geometry_model<3>, \
                                name, description); \
  }
  }
}

#endif
