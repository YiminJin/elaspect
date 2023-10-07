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


#include <elaspect/mesh_deformation/free_surface/projection.h>
#include <elaspect/mesh_deformation/handler.h>
#include <elaspect/gravity_model/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>

namespace elaspect
{
  namespace MeshDeformation
  {
    namespace FreeSurface
    {
      template <int dim>
      void 
      Projection<dim>::
      project_displacement_increment_onto_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                   const IndexSet &mesh_locally_owned,
                                                   const IndexSet &mesh_locally_relevant,
                                                   TrilinosWrappers::MPI::Vector &output) const
      {
        // stuff for iterating over the mesh
        const QGauss<dim-1> face_quadrature(mesh_deformation_dof_handler.get_fe().degree + 1);
        UpdateFlags update_flags(update_values | 
                                 update_quadrature_points |
                                 update_normal_vectors |
                                 update_JxW_values);
        FEFaceValues<dim> fs_fe_face_values(this->get_mapping(),
                                            mesh_deformation_dof_handler.get_fe(),
                                            face_quadrature,
                                            update_flags);
        FEFaceValues<dim> fe_face_values(this->get_mapping(),
                                         this->get_fe(),
                                         face_quadrature,
                                         update_flags);

        const unsigned int n_face_q_points = fe_face_values.n_quadrature_points;
        const unsigned int dofs_per_cell = fs_fe_face_values.dofs_per_cell;

        // stuff for assembling system
        std::vector<types::global_dof_index> cell_dof_indices(dofs_per_cell);
        Vector<double> cell_vector(dofs_per_cell);
        FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

        std::vector<Tensor<1,dim>> displacement_increments(n_face_q_points);
        

        // set up constraints
        AffineConstraints<double> mass_matrix_constraints(mesh_locally_relevant);
        DoFTools::make_hanging_node_constraints(mesh_deformation_dof_handler,
                                                mass_matrix_constraints);
        mass_matrix_constraints.close();

        // set up the matrix
        TrilinosWrappers::SparseMatrix mass_matrix;
        TrilinosWrappers::SparsityPattern sp(mesh_locally_owned,
                                             mesh_locally_owned,
                                             mesh_locally_relevant,
                                             this->get_mpi_communicator());
        DoFTools::make_sparsity_pattern(mesh_deformation_dof_handler, 
                                        sp, mass_matrix_constraints,
                                        false);
        sp.compress();
        mass_matrix.reinit(sp);

        FEValuesExtractors::Vector extractor(0);

        // make distributed vectors
        TrilinosWrappers::MPI::Vector rhs, dist_solution;
        rhs.reinit(mesh_locally_owned, this->get_mpi_communicator());
        dist_solution.reinit(mesh_locally_owned, this->get_mpi_communicator());

        const std::set<types::boundary_id> &fs_boundary_ids = 
          this->get_mesh_deformation_handler().get_free_surface_boundary_indicators();

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

              fs_cell->get_dof_indices(cell_dof_indices);
              fs_fe_face_values.reinit(fs_cell, face_no);
              fe_face_values.reinit(cell, face_no);
              fe_face_values[this->introspection().extractors.displacement].get_function_values(
                this->get_solution(), displacement_increments);

              cell_vector = 0;
              cell_matrix = 0;
              for (unsigned int q = 0; q < n_face_q_points; ++q)
              {
                // Select the direction onto which to project the displacement increments
                Tensor<1,dim> direction;
                if (projection_direction == SurfaceProjection::normal)
                  direction = fs_fe_face_values.normal_vector(q);
                else if (projection_direction == SurfaceProjection::vertical)
                  direction = this->get_gravity_model().gravity_vector(fs_fe_face_values.quadrature_point(q));
                else
                  AssertThrow(false, ExcInternalError());

                direction *= (direction.norm() > 0.0 ? 1. / direction.norm() : 0.0);

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    cell_matrix(i, j) += ( fs_fe_face_values[extractor].value(j, q) *
                                           fs_fe_face_values[extractor].value(i, q) ) *
                                         fs_fe_face_values.JxW(q);

                  cell_vector(i) += ( fs_fe_face_values[extractor].value(i, q) * direction ) *
                                    ( displacement_increments[q] * direction ) *
                                    fs_fe_face_values.JxW(q);
                }
              }

              mass_matrix_constraints.distribute_local_to_global(cell_matrix, cell_vector,
                                                                 cell_dof_indices,
                                                                 mass_matrix, rhs,
                                                                 false);
            }

        rhs.compress(VectorOperation::add);
        mass_matrix.compress(VectorOperation::add);

        // Jacobi seems to be fine here.  Other preconditioners (ILU, IC) run into troubles
        // because the matrix is mostly empty, since we don't touch internal vertices.
        TrilinosWrappers::PreconditionJacobi preconditioner_mass;
        preconditioner_mass.initialize(mass_matrix);

        SolverControl solver_control(rhs.size(), 
                                     this->get_parameters().mechanical_system_solver_tolerance * rhs.l2_norm());
        TrilinosWrappers::SolverCG cg(solver_control);
        cg.solve(mass_matrix, dist_solution, rhs, preconditioner_mass);

        mass_matrix_constraints.distribute(dist_solution);
        output = dist_solution;
      }


      template <int dim>
      void
      Projection<dim>::
      make_boundary_constraints(const TrilinosWrappers::MPI::Vector &mesh_displacements,
                                const DoFHandler<dim>               &mesh_deformation_dof_handler,
                                AffineConstraints<double>           &mesh_deformation_constraints,
                                const std::set<types::boundary_id>  &fs_boundary_ids) const
      {
        // For the free surface indicators we constrain the displacement to be u.n
        const IndexSet mesh_locally_owned = mesh_deformation_dof_handler.locally_owned_dofs();
        IndexSet mesh_locally_relevant;
        DoFTools::extract_locally_relevant_dofs(mesh_deformation_dof_handler, 
                                                mesh_locally_relevant);
        
        TrilinosWrappers::MPI::Vector boundary_displacement_increment(
          mesh_locally_owned, mesh_locally_relevant, this->get_mpi_communicator());

        project_displacement_increment_onto_boundary(mesh_deformation_dof_handler, 
                                                     mesh_locally_owned,
                                                     mesh_locally_relevant,
                                                     boundary_displacement_increment);

        // insert the relevant part of the solution into the mesh constraints
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


      template <int dim>
      void Projection<dim>::declare_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Mesh deformation");
        {
          prm.enter_subsection("Free surface");
          {
            prm.enter_subsection("Projection");
            {
              prm.declare_entry("Projection direction", "normal",
                                Patterns::Selection("normal|vertical"),
                                "");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void Projection<dim>::parse_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Mesh deformation");
        {
          prm.enter_subsection("Free surface");
          {
            prm.enter_subsection("Projection");
            {
              std::string projection_dir = prm.get("Projection direction");

              if (projection_dir == "normal")
                projection_direction = SurfaceProjection::normal;
              else if (projection_dir == "vertical")
                projection_direction = SurfaceProjection::vertical;
              else
                AssertThrow(false, ExcNotImplemented());
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
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
      ELASPECT_REGISTER_FREE_SURFACE_MODEL(Projection,
                                       "projection",
                                       "")
    }
  }
}
