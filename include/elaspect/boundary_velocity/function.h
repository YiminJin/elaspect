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


#ifndef _elaspect_boundary_velocity_function_h
#define _elaspect_boundary_velocity_function_h

#include <elaspect/boundary_velocity/interface.h>
#include <elaspect/simulator_access.h>
#include <elaspect/utilities.h>

#include <deal.II/base/parsed_function.h>

namespace elaspect
{
  namespace BoundaryVelocity
  {
    using namespace dealii;

    /**
     * A class that implements velocity boundary conditions based on a
     * functional description provided in the input file.
     *
     * @ingroup BoundaryVelocities
     */
    template <int dim>
    class Function : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        Function ();

        /**
         * Return the boundary velocity as a function of position. For the
         * current class, this function obviously simply returns a zero
         * tensor.
         */
        Tensor<1,dim>
        boundary_velocity (const types::boundary_id boundary_indicator,
                           const Point<dim> &position) const override;

        // avoid -Woverloaded-virtual warning until the deprecated function
        // is removed from the interface:
        using Interface<dim>::boundary_velocity;

        /**
         * A function that is called at the beginning of each time step to
         * indicate what the model time is for which the boundary values will
         * next be evaluated. For the current class, the function passes to
         * the parsed function what the current time is.
         */
        void
        update () override;

        /**
         * Declare the parameters this class takes through input files. The
         * default implementation of this function does not describe any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * The default implementation of this function does not read any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * A function object representing the components of the velocity.
         */
        Functions::ParsedFunction<dim> boundary_velocity_function;

        /**
         * The coordinate representation to evaluate the function. Possible
         * choices are depth, cartesian and spherical.
         */
        Utilities::Coordinates::CoordinateSystem coordinate_system;

        /**
         * Whether to specify velocity in x, y, z components, or
         * r, phi, theta components.
         */
        bool use_spherical_unit_vectors;
    };
  }
}


#endif
