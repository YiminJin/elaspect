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


#include <elaspect/initial_topography/function.h>
#include <elaspect/geometry_model/box.h>
#include <elaspect/geometry_model/sphere.h>

namespace elaspect
{
  namespace InitialTopography
  {
    template <int dim>
    Function<dim>::Function()
      :
      max_topo(0),
      initial_topography_function(1),
      coordinate_system(Utilities::Coordinates::CoordinateSystem::cartesian)
    {}


    template <int dim>
    double
    Function<dim>::
    value(const Point<dim-1> &surface_point) const
    {
      Point<dim> global_point;
      if (Plugins::plugin_type_matches<GeometryModel::Box<dim>>(this->get_geometry_model()))
      {
        for (unsigned int d = 0; d < dim-1; ++d)
          global_point[d] = surface_point[d];
      }
      else if (Plugins::plugin_type_matches<GeometryModel::Sphere<dim>>(this->get_geometry_model()))
      {
        // TODO set radius
      }
      else
        AssertThrow(false, ExcNotImplemented());

      const double topo = initial_topography_function.value(global_point);

      return topo;
    }


    template <int dim>
    double
    Function<dim>::max_topography() const
    {
      return max_topo;
    }


    template <int dim>
    void
    Function<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Initial topography");
      {
        prm.enter_subsection("Function");
        {
          prm.declare_entry("Maximum topography value", "0",
                            Patterns::Double(0),
                            "The maximum value the topography given by "
                            "the function can take.");
          prm.declare_entry("Coordinate system", "cartesian",
                            Patterns::Selection("cartesian|spherical"),
                            "A selection that determines the assumed coordinate "
                            "system for the function variables. Allowed values "
                            "are `cartesian' an `spherical'. `spherical' coordinates "
                            "are interpreted as r,phi or r,phi,theta in 2d/3d "
                            "respectively with theta being the polar angle.");

          Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Function<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Initial topography");
      {
        prm.enter_subsection("Function");
        {
          coordinate_system = Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));
          max_topo = prm.get_double("Maximum topography value");

          try
          {
            initial_topography_function.parse_parameters(prm);
          }
          catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Boundary traction model.Function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'"
                      << "More information about the cause of the parse error \n"
                      << "is shown below.\n";
            throw;
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
  namespace InitialTopography
  {
    ELASPECT_REGISTER_INITIAL_TOPOGRAPHY_MODEL(Function,
                                           "function",
                                           "")
  }
}
