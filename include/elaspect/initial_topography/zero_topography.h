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


#ifndef _elaspect_initial_topography_zero_topography_h
#define _elaspect_initial_topography_zero_topography_h

#include <elaspect/initial_topography/interface.h>

namespace elaspect
{
  namespace InitialTopography
  {
    template <int dim>
    class ZeroTopography : public Interface<dim>
    {
      public:
        double
        value (const Point<dim-1> &p) const override;

        double max_topography () const override;
    };
  }
}

#endif
