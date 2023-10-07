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


#ifndef _elaspect_boundary_temperature_function_h
#define _elaspect_boundary_temperature_function_h

#include <elaspect/boundary_temperature/interface.h>
#include <elaspect/simulator_access.h>
#include <elaspect/utilities.h>

#include <deal.II/base/parsed_function.h>

namespace elaspect
{
  namespace BoundaryTemperature
  {
    using namespace dealii;

    template <int dim>
    class Function : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        Function ();

        void update () override;

        double
        boundary_temperature (const types::boundary_id boundary_indicator,
                              const Point<dim> &position) const override;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        Functions::ParsedFunction<dim> boundary_temperature_function;

        Utilities::Coordinates::CoordinateSystem coordinate_system;
    };
  }
}

#endif
