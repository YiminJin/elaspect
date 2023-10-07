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


#include <elaspect/simulator.h>

namespace elaspect
{
  template <int dim>
  void
  Simulator<dim>::assemble_and_solve_thermo_system ()
  {
    if (parameters.include_heat_transport)
    {
      assemble_thermo_system();
      solve_thermo_system();
    }
  }


  template <int dim>
  void 
  Simulator<dim>::assemble_and_solve_hydro_system ()
  {

  }


  template <int dim>
  bool
  Simulator<dim>::assemble_and_solve_mechanical_system ()
  {
    const unsigned int block_idx = introspection.block_indices.displacement;

    TrilinosWrappers::MPI::BlockVector incremental_displacement(
      introspection.index_sets.system_partitioning, mpi_communicator);

    bool converged = false;
    if (parameters.constitutive_relation & ConstitutiveRelation::plasticity)
    {
      TrilinosWrappers::MPI::BlockVector 
      newton_update(introspection.index_sets.system_partitioning, 
                    mpi_communicator),
      incremental_displacement_backup(introspection.index_sets.system_partitioning,
                                      introspection.index_sets.system_relevant_partitioning, 
                                      mpi_communicator);

      double initial_residual = 0;
      nonlinear_iteration = 0;
      while (true)
      {
        if (nonlinear_iteration == 1)
          compute_current_constraints();

        assemble_mechanical_system();

        double residual = system_rhs.block(block_idx).l2_norm();
        if (nonlinear_iteration == 0)
          initial_residual = residual;

        pcout << "   iteration " << nonlinear_iteration
              << ": residual = " << residual
              << std::endl;

        // break the loop if the nonlinear system is converged...
        if (residual / initial_residual < parameters.nonlinear_tolerance)
        {
          converged = true;
          break;
        }

        // or the number of nonlinear iterations has reached the upper limit.
        if (nonlinear_iteration >= parameters.max_nonlinear_iterations)
        {
          converged = false;
          break;
        }

        solve_mechanical_system(newton_update);

        if (nonlinear_iteration == 0)
        {
          solution.block(block_idx) = newton_update.block(block_idx);
          apply_return_mapping();
        }
        else
        {
          // execute line search
          incremental_displacement_backup.block(block_idx) = solution.block(block_idx);

          double alpha = 1;
          unsigned int line_search_step = 0;

          reassemble_tangent_matrix = false;
          while (line_search_step < parameters.max_line_search_steps)
          {
            incremental_displacement.block(block_idx) = incremental_displacement_backup.block(block_idx);
            incremental_displacement.block(block_idx).add(alpha, newton_update.block(block_idx));
            solution.block(block_idx) = incremental_displacement.block(block_idx);

            apply_return_mapping();
            assemble_mechanical_system();
            double test_residual = system_rhs.block(block_idx).l2_norm();

            if (test_residual * test_residual <=
                (1 - 2 * parameters.line_search_beta * alpha) * residual * residual)
              break;
            else
              alpha *= parameters.line_search_delta;

            ++line_search_step;
          }

          pcout << "   line search: " << line_search_step
                << " steps." << std::endl;

          reassemble_tangent_matrix = true;
        }

        if (parameters.run_postprocessors_on_nonlinear_iterations)
          postprocess();

        ++nonlinear_iteration;
      }
    }
    else
    {
      assemble_mechanical_system();
      solve_mechanical_system(incremental_displacement);
      solution.block(block_idx) = incremental_displacement.block(block_idx);

      converged = true;
    }

    return converged;
  }


  template <int dim>
  void
  Simulator<dim>::assemble_and_solve_qpd_system ()
  {
    TimerOutput::Scope timer(computing_timer, 
                             (parameters.use_ALE_method ? "QPD advection" : "QPD projection"));

    pcout << "   "
          << (parameters.use_ALE_method ? "Advecting" : "Projecting")
          << " quadrature point data " 
          << (parameters.use_ALE_method ? "across" : "onto")
          << " the grid... "
          << std::flush;

    // The QPD fields can be divided into two categories: compositional fields and
    // physical fields, which are discretized with different finite element spaces.
    std::vector<QPDField> compositional_fields, physical_fields;

    for (unsigned int i = 0; i < introspection.n_compositional_fields; ++i)
      compositional_fields.push_back(QPDField(introspection, "composition", i));

    // Physical fields that require projection/advection include stress components
    // and plastic strain.
    for (unsigned int i = 0; i < SymmetricTensor<2,dim>::n_independent_components; ++i)
      physical_fields.push_back(QPDField(introspection, "stress", i));

    if (parameters.constitutive_relation & ConstitutiveRelation::plasticity)
      physical_fields.push_back(QPDField(introspection, "plastic strain", 0));

    // Assemble and solve the two categories of fields seperatedly.
    assemble_qpd_system(compositional_fields);
    solve_qpd_system(compositional_fields);

    assemble_qpd_system(physical_fields);
    solve_qpd_system(physical_fields);

    pcout << "done." << std::endl;
  }
}


// explicit instantiations
namespace elaspect
{
#define INSTANTIATE(dim) \
  template void Simulator<dim>::assemble_and_solve_thermo_system(); \
  template void Simulator<dim>::assemble_and_solve_hydro_system(); \
  template bool Simulator<dim>::assemble_and_solve_mechanical_system(); \
  template void Simulator<dim>::assemble_and_solve_qpd_system();

  ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
