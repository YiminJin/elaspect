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


#ifndef _elaspect_mesh_deformation_mesh_smoothing_laplacian_h
#define _elaspect_mesh_deformation_mesh_smoothing_laplacian_h

#include <elaspect/mesh_deformation/mesh_smoothing/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace MeshDeformation
  {
    namespace MeshSmoothing
    {
      template <int dim>
      class Laplacian : public Interface<dim>, public SimulatorAccess<dim>
      {
        public:
          void
          compute_mesh_displacements(const DoFHandler<dim>           &dof_handler,
                                     const AffineConstraints<double> &constraints,
                                     TrilinosWrappers::MPI::Vector   &displacements) const override;
      };
    }
  }
}

#endif
