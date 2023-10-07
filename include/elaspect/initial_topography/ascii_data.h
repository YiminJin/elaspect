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


#ifndef _elaspect_initial_topography_ascii_data_h
#define _elaspect_initial_topography_ascii_data_h

#include <elaspect/initial_topography/interface.h>
#include <elaspect/structured_data.h>

namespace elaspect
{
  namespace InitialTopography
  {
    template <int dim>
    class AsciiData : public Utilities::AsciiDataBoundary<dim>, public Interface<dim>
    {
      public:
        AsciiData();

        void
        initialize() override;

        using Utilities::AsciiDataBoundary<dim>::initialize;

        double
        value(const Point<dim-1> &surface_point) const override;

        double max_topography() const override;

        static
        void
        declare_parameters(ParameterHandler &prm);

        void
        parse_parameters(ParameterHandler &prm) override;

      private:
        types::boundary_id surface_boundary_id;
    };
  }
}

#endif
