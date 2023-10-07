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


#include <elaspect/mesh_deformation/free_surface/simple.h>
#include <elaspect/mesh_deformation/handler.h>

#include <deal.II/dofs/dof_tools.h>

namespace elaspect
{
  namespace MeshDeformation
  {
    namespace FreeSurface
    {
      template <int dim>
      void
      Simple<dim>::
      make_boundary_constraints(const TrilinosWrappers::MPI::Vector &mesh_displacements,
                                const DoFHandler<dim>               &mesh_deformation_dof_handler,
                                AffineConstraints<double>           &mesh_deformation_constraints,
                                const std::set<types::boundary_id>  &fs_boundary_ids) const
      {
        const IndexSet mesh_locally_owned = mesh_deformation_dof_handler.locally_owned_dofs();
        IndexSet mesh_locally_relevant;
        DoFTools::extract_locally_relevant_dofs(mesh_deformation_dof_handler,
                                                mesh_locally_relevant);
        
        TrilinosWrappers::MPI::Vector boundary_displacement_increment, distributed_vector;
        boundary_displacement_increment.reinit(mesh_locally_owned, 
                                               mesh_locally_relevant, 
                                               this->get_mpi_communicator());
        distributed_vector.reinit(mesh_locally_owned, 
                                  this->get_mpi_communicator());

        std::vector<types::global_dof_index> vertex_dof_indices(dim), mesh_vertex_dof_indices(dim);
        std::vector<double> vertex_dof_values(dim);

        std::vector<bool> vertex_visited(this->get_triangulation().n_vertices(), false);

        typename DoFHandler<dim>::active_cell_iterator
        cell = this->get_dof_handler().begin_active(),
        endc = this->get_dof_handler().end();
        typename DoFHandler<dim>::active_cell_iterator
        fs_cell = mesh_deformation_dof_handler.begin_active();
        for (; cell != endc; ++cell, ++fs_cell)
          if (cell->at_boundary() && cell->is_locally_owned())
            for (const unsigned int face_no : cell->face_indices())
            {
              const types::boundary_id boundary_id = cell->face(face_no)->boundary_id();
              if (fs_boundary_ids.find(boundary_id) == fs_boundary_ids.end())
                continue;

              for (const unsigned int vertex_no : cell->vertex_indices())
              {
                if (vertex_visited[cell->vertex_index(vertex_no)])
                  continue;

                for (unsigned int d = 0; d < dim; ++d)
                {
                  vertex_dof_indices[d] = cell->vertex_dof_index(vertex_no, d);
                  mesh_vertex_dof_indices[d] = fs_cell->vertex_dof_index(vertex_no, d);
                }

                this->get_solution().extract_subvector_to(vertex_dof_indices, vertex_dof_values);
                distributed_vector.set(mesh_vertex_dof_indices, vertex_dof_values);

                vertex_visited[cell->vertex_index(vertex_no)] = true;
              }
            }

        distributed_vector.compress(VectorOperation::insert);
        boundary_displacement_increment = distributed_vector;

        const IndexSet constrained_dofs = 
          DoFTools::extract_boundary_dofs(mesh_deformation_dof_handler,
                                          ComponentMask(dim, true),
                                          fs_boundary_ids);

        for (unsigned int i = 0; i < constrained_dofs.n_elements(); ++i)
        {
          types::global_dof_index index = constrained_dofs.nth_index_in_set(i);
          if (mesh_deformation_constraints.can_store_line(index))
            if (mesh_deformation_constraints.is_constrained(index) == false)
            {
              mesh_deformation_constraints.add_line(index);
              mesh_deformation_constraints.set_inhomogeneity(
                index, boundary_displacement_increment[index] + mesh_displacements[index]);
            }
        }
      }
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace MeshDeformation
  {
    namespace FreeSurface
    {
      ELASPECT_REGISTER_FREE_SURFACE_MODEL(Simple,
                                       "simple",
                                       "")
    }
  }
}
