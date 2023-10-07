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


#include <elaspect/gravity_model/radial.h>
#include <elaspect/geometry_model/interface.h>

namespace elaspect
{
  namespace GravityModel
  {
// ------------------------------ RadialConstant -------------------

    template <int dim>
    Tensor<1,dim>
    RadialConstant<dim>::gravity_vector (const Point<dim> &p) const
    {
      if (p.norm() == 0.0)
        return Tensor<1,dim>();

      const double r = p.norm();
      return -magnitude * p/r;
    }


    template <int dim>
    void
    RadialConstant<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Radial constant");
        {
          prm.declare_entry ("Magnitude", "9.81",
                             Patterns::Double (),
                             "Magnitude of the gravity vector in $m/s^2$. For positive values "
                             "the direction is radially inward towards the center of the earth.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    RadialConstant<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Radial constant");
        {
          magnitude = prm.get_double ("Magnitude");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

      AssertThrow (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::spherical ||
                   this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::ellipsoidal,
                   ExcMessage ("Gravity model 'radial constant' should not be used with geometry models that "
                               "do not have either a spherical or ellipsoidal natural coordinate system."));
    }


// ----------------------------- RadialLinear ----------------------

    template <int dim>
    Tensor<1,dim>
    RadialLinear<dim>::gravity_vector (const Point<dim> &p) const
    {
      if (p.norm() == 0.0)
        return Tensor<1,dim>();

      const double depth = this->get_geometry_model().depth(p);
      const double maximal_depth = this->get_geometry_model().maximal_depth();

      return  (-magnitude_at_surface * (1.0 - depth/maximal_depth)
               -magnitude_at_bottom  * (depth/maximal_depth))
              * p/p.norm();
    }


    template <int dim>
    void
    RadialLinear<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Radial linear");
        {
          prm.declare_entry ("Magnitude at surface", "9.8",
                             Patterns::Double (),
                             "Magnitude of the radial gravity vector "
                             "at the surface of the domain. "
                             "Units: \\si{\\meter\\per\\second\\squared}.");
          prm.declare_entry ("Magnitude at bottom", "10.7",
                             Patterns::Double (),
                             "Magnitude of the radial gravity vector "
                             "at the bottom of the domain. `Bottom' means the"
                             "maximum depth in the chosen geometry, and for "
                             "example represents the core-mantle boundary in "
                             "the case of the `spherical shell' geometry model, "
                             "and the center in the case of the `sphere' "
                             "geometry model. Units: \\si{\\meter\\per\\second\\squared}.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    RadialLinear<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Radial linear");
        {
          magnitude_at_surface = prm.get_double ("Magnitude at surface");
          magnitude_at_bottom = prm.get_double ("Magnitude at bottom");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

      AssertThrow (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::spherical ||
                   this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::ellipsoidal,
                   ExcMessage ("Gravity model 'radial linear' should not be used with geometry models that "
                               "do not have either a spherical or ellipsoidal natural coordinate system."));

    }


// ------------------------------ RadialConstant -------------------

    template <int dim>
    Tensor<1,dim>
    RadialFunction<dim>::gravity_vector (const Point<dim> &p) const
    {
      if (p.norm() == 0.0)
        return Tensor<1,dim>();

      const double r = p.norm();
      return -magnitude_function.value(Point<1>(r)) * p / r;
    }


    template <int dim>
    void
    RadialFunction<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Radial function");
        {
          Functions::ParsedFunction<1>::declare_parameters (prm, 1);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    RadialFunction<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Radial function");
        {
          try
          {
            magnitude_function.parse_parameters (prm);
          }
          catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Gravity model.Radial function'\n"
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
  namespace GravityModel
  {
    ELASPECT_REGISTER_GRAVITY_MODEL(RadialConstant,
                                "radial constant",
                                "A gravity model in which the gravity has a constant magnitude "
                                "and the direction is radial (pointing inward if the value "
                                "is positive). The magnitude is read from the parameter "
                                "file in subsection 'Radial constant'.")

    ELASPECT_REGISTER_GRAVITY_MODEL(RadialLinear,
                                "radial linear",
                                "A gravity model which is radial (pointing inward if the gravity "
                                "is positive) and the magnitude changes linearly with depth. The "
                                "magnitude of gravity at the surface and bottom is read from the "
                                "input file in a section ``Gravity model/Radial linear''.")

    ELASPECT_REGISTER_GRAVITY_MODEL(RadialFunction,
                                "radial function",
                                "A gravity model which is radial (pointing inward if the gravity "
                                "is positive) and the magnitude is described by a function.")
  }
}
