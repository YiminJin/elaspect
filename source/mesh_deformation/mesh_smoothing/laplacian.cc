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


#include <elaspect/mesh_deformation/mesh_smoothing/laplacian.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

namespace elaspect
{
  namespace MeshDeformation
  {
    namespace MeshSmoothing
    {
      template <int dim>
      void
      Laplacian<dim>::
      compute_mesh_displacements(const DoFHandler<dim>           &dof_handler,
                                 const AffineConstraints<double> &constraints,
                                 TrilinosWrappers::MPI::Vector   &displacements) const
      {
        const FiniteElement<dim> &fe = dof_handler.get_fe();
        QGauss<dim> quadrature(fe.degree + 1);
        UpdateFlags update_flags(update_values | update_JxW_values | update_gradients);
        FEValues<dim> fe_values(fe, quadrature, update_flags);

        const unsigned int dofs_per_cell = fe_values.dofs_per_cell;
        const unsigned int n_q_points = fe_values.n_quadrature_points;

        std::vector<types::global_dof_index> cell_dof_indices(dofs_per_cell);
        FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
        Vector<double> cell_rhs(dofs_per_cell);

        // We are just solving a Laplacian in each spatial direction, so
        // the degrees of freedom for different dimensions do not couple.
        Table<2,DoFTools::Coupling> coupling(dim, dim);
        coupling.fill(DoFTools::none);

        for (unsigned int c = 0; c < dim; ++c)
          coupling[c][c] = DoFTools::always;

        IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
        IndexSet locally_relevant_dofs;
        DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

        TrilinosWrappers::SparsityPattern sp(locally_owned_dofs,
                                             locally_owned_dofs,
                                             locally_relevant_dofs,
                                             this->get_mpi_communicator());

        DoFTools::make_sparsity_pattern(dof_handler, coupling, sp,
                                        constraints, false,
                                        Utilities::MPI::
                                        this_mpi_process(this->get_mpi_communicator()));
        sp.compress();

        TrilinosWrappers::SparseMatrix system_matrix;
        system_matrix.reinit(sp);

        // carry out the solution
        FEValuesExtractors::Vector extractor(0);

        TrilinosWrappers::MPI::Vector system_rhs(locally_owned_dofs, this->get_mpi_communicator());
        TrilinosWrappers::MPI::Vector solution(locally_owned_dofs, this->get_mpi_communicator());

        for (const auto &cell : dof_handler.active_cell_iterators())
          if (cell->is_locally_owned())
          {
            cell->get_dof_indices(cell_dof_indices);
            fe_values.reinit(cell);

            cell_matrix = 0;
            cell_rhs = 0;
            for (unsigned int q = 0; q < n_q_points; ++q)
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  cell_matrix(i, j) += scalar_product(fe_values[extractor].gradient(i, q),
                                                      fe_values[extractor].gradient(j, q))
                                       * fe_values.JxW(q);
            
            constraints.distribute_local_to_global(cell_matrix, cell_rhs,
                                                   cell_dof_indices, 
                                                   system_matrix, system_rhs,
                                                   false);
          }

        system_rhs.compress(VectorOperation::add);
        system_matrix.compress(VectorOperation::add);

        const double tolerance = system_rhs.l2_norm() * 1e-8;
        if (tolerance < 1e-50)
          return;

        this->get_pcout() << "   Executing Laplacian smoothing... " << std::flush;

        // Make the AMG preconditioner
        std::vector<std::vector<bool>> constant_modes;
        DoFTools::extract_constant_modes(dof_handler, 
                                          ComponentMask(dim, true),
                                          constant_modes);

        TrilinosWrappers::PreconditionAMG preconditioner;
        TrilinosWrappers::PreconditionAMG::AdditionalData Amg_data;

        Amg_data.constant_modes = constant_modes;
        Amg_data.elliptic = true;
        Amg_data.higher_order_elements = false;
        Amg_data.smoother_sweeps = 2;
        Amg_data.aggregation_threshold = 0.02;

        preconditioner.initialize(system_matrix);

        SolverControl solver_control(system_rhs.size(), 1e-8 * system_rhs.l2_norm());
        TrilinosWrappers::SolverCG cg(solver_control);

        cg.solve(system_matrix, solution, system_rhs, preconditioner);
        
        this->get_pcout() << "done." << std::endl;

        constraints.distribute(solution);
        displacements = solution;
      }
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace MeshDeformation
  {
    namespace MeshSmoothing
    {
      ELASPECT_REGISTER_MESH_SMOOTHING_MODEL(Laplacian,
                                         "laplacian",
                                         "")
    }
  }
}
