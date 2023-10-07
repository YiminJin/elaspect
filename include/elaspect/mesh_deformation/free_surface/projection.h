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


#ifndef _elaspect_mesh_deformation_free_surface_projection_h
#define _elaspect_mesh_deformation_free_surface_projection_h

#include <elaspect/mesh_deformation/free_surface/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace MeshDeformation
  {
    namespace FreeSurface
    {
      template <int dim>
      class Projection : public Interface<dim>, public SimulatorAccess<dim>
      {
        public:
          void 
          make_boundary_constraints(const TrilinosWrappers::MPI::Vector &mesh_displacements,
                                    const DoFHandler<dim>               &mesh_deformation_dof_handler,
                                    AffineConstraints<double>           &mesh_deformation_constraints,
                                    const std::set<types::boundary_id>  &boundary_id) const override;

          static
          void declare_parameters(ParameterHandler &prm);

          void parse_parameters(ParameterHandler &prm) override;

        private:
          void 
          project_displacement_increment_onto_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                       const IndexSet &mesh_locally_owned,
                                                       const IndexSet &mesh_locally_relevant,
                                                       TrilinosWrappers::MPI::Vector &output) const;

          struct SurfaceProjection
          {
            enum Direction { normal, vertical };
          };

          typename SurfaceProjection::Direction projection_direction;
      };
    }
  }
}

#endif
