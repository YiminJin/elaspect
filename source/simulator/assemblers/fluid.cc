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


#include <elaspect/simulator/assemblers/hydro_system.h>

namespace elaspect
{
  namespace Assemblers
  {
    template <int dim>
    void
    HydroUnsaturatedFlow<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &/*scratch_base*/,
             internal::Assembly::CopyData::CopyDataBase<dim> &/*data_base*/) const
    {

    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace Assemblers
  {
#define INSTANTIATE(dim) \
    template class HydroUnsaturatedFlow<dim>;

    ELASPECT_INSTANTIATE(INSTANTIATE)
  }
}
