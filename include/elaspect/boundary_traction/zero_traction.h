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


#ifndef _elaspect_boundary_traction_zero_traction_h
#define _elaspect_boundary_traction_zero_traction_h

#include <elaspect/boundary_traction/interface.h>

namespace elaspect
{
  namespace BoundaryTraction
  {
    using namespace dealii;

    template <int dim>
    class ZeroTraction : public Interface<dim>
    {
      public:
        Tensor<1,dim>
        boundary_traction (const types::boundary_id boundary_indicator,
                           const Point<dim> &position,
                           const Tensor<1,dim> &normal_vector) const override;
    };
  }
}

#endif
