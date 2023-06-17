#include <elaspect/mesh_deformation/mesh_smoothing/hyperelasticity.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/physics/elasticity/kinematics.h>
#include <deal.II/physics/elasticity/standard_tensors.h>

namespace elaspect
{
  namespace MeshDeformation
  {
    namespace MeshSmoothing
    {
      namespace internal
      {
        template <int dim>
        SymmetricTensor<4,dim> dC_inv_dC (const SymmetricTensor<2,dim> &C_inv)
        {
          SymmetricTensor<4,dim> result;
          for (unsigned int A=0; A<dim; ++A)
            for (unsigned int B=A; B<dim; ++B)
              for (unsigned int C=0; C<dim; ++C)
                for (unsigned int D=C; D<dim; ++D)
                  result[A][B][C][D] -= 0.5 * (C_inv[A][C] * C_inv[B][D] + C_inv[A][D] * C_inv[B][C]);

          return result;
        }
      }


      template <int dim>
      void
      Hyperelasticity<dim>::
      make_newton_constraints(const AffineConstraints<double>     &constraints,
                              AffineConstraints<double>           &newton_constraints,
                              const TrilinosWrappers::MPI::Vector &displacements,
                              const bool                           use_inhomogeneity) const
      {
        const IndexSet local_lines = constraints.get_local_lines();
        newton_constraints.clear();
        newton_constraints.reinit(local_lines);

        newton_constraints.merge(constraints);

        if (use_inhomogeneity)
        {
          for (auto index = local_lines.begin(); index != local_lines.end(); ++index)
            if (newton_constraints.is_inhomogeneously_constrained(*index))
              newton_constraints.set_inhomogeneity(
                  *index, newton_constraints.get_inhomogeneity(*index) - displacements[*index]);
        }
        else
        {
          for (auto index = local_lines.begin(); index != local_lines.end(); ++index)
            if (newton_constraints.is_constrained(*index))
              newton_constraints.set_inhomogeneity(*index, 0);
        }

        newton_constraints.close();
      }


      template <int dim>
      void
      Hyperelasticity<dim>::
      assemble_linear_system(const DoFHandler<dim>               &dof_handler,
                             const AffineConstraints<double>     &newton_constraints,
                             const TrilinosWrappers::MPI::Vector &displacements,
                             TrilinosWrappers::SparseMatrix      &system_matrix,
                             TrilinosWrappers::MPI::Vector       &system_rhs,
                             const bool                           reassemble_matrix) const
      {
        const FiniteElement<dim> &fe = dof_handler.get_fe();
        QGauss<dim> quadrature(fe.degree + 1);
        UpdateFlags update_flags(update_gradients | update_JxW_values);
        FEValues<dim> fe_values(fe, quadrature, update_flags);

        const unsigned int n_q_points = fe_values.n_quadrature_points;
        const unsigned int dofs_per_cell = fe_values.dofs_per_cell;

        const FEValuesExtractors::Vector extractor(0);

        std::vector<types::global_dof_index> cell_dof_indices(dofs_per_cell);
        FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
        Vector<double> cell_rhs(dofs_per_cell);

        std::vector<Tensor<2,dim>>          Grad_phi(dofs_per_cell);
        std::vector<SymmetricTensor<2,dim>> symm_FT_Grad_phi(dofs_per_cell);
        std::vector<Tensor<2,dim>>          u_gradients(n_q_points);

        double ln_J;
        Tensor<2,dim> F;
        SymmetricTensor<2,dim> C, C_inv, S;
        SymmetricTensor<4,dim> D;

        system_rhs = 0;
        if (reassemble_matrix)
          system_matrix = 0;

        for (const auto &cell : dof_handler.active_cell_iterators())
          if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);
            fe_values[extractor].get_function_gradients(displacements, u_gradients);

            cell_rhs = 0;
            if (reassemble_matrix)
              cell_matrix = 0;
            cell->get_dof_indices(cell_dof_indices);

            for (unsigned int q = 0; q < n_q_points; ++q)
            {
              F     = Physics::Elasticity::Kinematics::F(u_gradients[q]);
              C     = Physics::Elasticity::Kinematics::C(F);
              C_inv = invert(C);
              ln_J  = std::log(determinant(F));
              S     = (Physics::Elasticity::StandardTensors<dim>::I - C_inv)
                      + (lambda * ln_J) * C_inv;
              D     = lambda * outer_product(C_inv, C_inv) -
                      (2.0 * (1.0 - lambda * ln_J)) * internal::dC_inv_dC(C_inv);

              for (unsigned int k = 0; k < dofs_per_cell; ++k)
              {
                Grad_phi[k] = fe_values[extractor].gradient(k, q);
                symm_FT_Grad_phi[k] = symmetrize(transpose(F) * Grad_phi[k]);               
              }

              const double JxW = fe_values.JxW(q);

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                cell_rhs(i) -= (S * symm_FT_Grad_phi[i]) * JxW;

                if (reassemble_matrix)
                {
                  // Since the cell matrix is symmetric, we can exploit this property by
                  // building only the lower half of the cell matrix and copying the
                  // values to the upper half.
                  const unsigned int comp_i = fe.system_to_component_index(i).first;
                  for (unsigned int j = 0; j <= i; ++j)
                  {
                    const unsigned int comp_j = fe.system_to_component_index(j).first;
                    // The material contribution:
                    cell_matrix(i, j) += (symm_FT_Grad_phi[i] * D * symm_FT_Grad_phi[j]) * JxW;
                    // The geometrical contribution:
                    if (comp_i == comp_j)
                      cell_matrix(i, j) += (Grad_phi[i][comp_i] * S * Grad_phi[j][comp_j]) * JxW;
                  }
                }
              }
            }

            // Copy the lower half of the cell matrix into the upper half:
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              for (unsigned int j = i + 1; j < dofs_per_cell; ++j)
                cell_matrix(i, j) = cell_matrix(j, i);

            newton_constraints.distribute_local_to_global(cell_matrix,
                                                          cell_rhs,
                                                          cell_dof_indices,
                                                          system_matrix,
                                                          system_rhs);
          }

        system_rhs.compress(VectorOperation::add);
        if (reassemble_matrix)
          system_matrix.compress(VectorOperation::add);
      }


      template <int dim>
      void
      Hyperelasticity<dim>::
      solve_linear_system(const DoFHandler<dim>                &dof_handler,
                          const AffineConstraints<double>      &newton_constraints,
                          const TrilinosWrappers::SparseMatrix &system_matrix,
                          const TrilinosWrappers::MPI::Vector  &system_rhs,
                          TrilinosWrappers::MPI::Vector        &newton_update) const
      {
        // Make the AMG preconditioner
        std::vector<std::vector<bool>> constant_modes;
        DoFTools::extract_constant_modes(dof_handler,
                                         ComponentMask(dim, true),
                                         constant_modes);

        TrilinosWrappers::PreconditionAMG::AdditionalData Amg_data;
        Amg_data.constant_modes        = constant_modes;
        Amg_data.elliptic              = true;
        Amg_data.higher_order_elements = false;
        Amg_data.smoother_sweeps       = 2;
        Amg_data.aggregation_threshold = 0.02;

        TrilinosWrappers::PreconditionAMG preconditioner;
        preconditioner.initialize(system_matrix);

        SolverControl solver_control(system_rhs.size(),
                                     linear_solver_tolerance * system_rhs.l2_norm());
        TrilinosWrappers::SolverCG cg(solver_control);

        newton_constraints.set_zero(newton_update);

        cg.solve(system_matrix, newton_update, system_rhs, preconditioner);

        newton_constraints.distribute(newton_update);
      }


      template <int dim>
      void
      Hyperelasticity<dim>::
      compute_mesh_displacements(const DoFHandler<dim>           &dof_handler,
                                 const AffineConstraints<double> &constraints,
                                 TrilinosWrappers::MPI::Vector   &displacements) const
      {
        this->get_pcout() << "   Executing hyperelastic smoothing... " << std::flush;

        // set up the system tangent matrix
        const IndexSet &locally_owned_dofs = dof_handler.locally_owned_dofs();
        const IndexSet &locally_relevant_dofs = constraints.get_local_lines();

        TrilinosWrappers::SparsityPattern sp(locally_owned_dofs,
                                             locally_owned_dofs,
                                             locally_relevant_dofs,
                                             this->get_mpi_communicator());

        DoFTools::make_sparsity_pattern(dof_handler, sp,
                                        constraints, false,
                                        Utilities::MPI::
                                        this_mpi_process(this->get_mpi_communicator()));
        sp.compress();

        TrilinosWrappers::SparseMatrix system_matrix;
        system_matrix.reinit(sp);

        // set up vectors
        TrilinosWrappers::MPI::Vector system_rhs(locally_owned_dofs, this->get_mpi_communicator());
        TrilinosWrappers::MPI::Vector newton_update(locally_owned_dofs, this->get_mpi_communicator());
        TrilinosWrappers::MPI::Vector test_displacements(locally_owned_dofs, locally_relevant_dofs, this->get_mpi_communicator());
        TrilinosWrappers::MPI::Vector dist_displacements(locally_owned_dofs, this->get_mpi_communicator());

        // set up constraints
        AffineConstraints<double> newton_constraints;
        newton_constraints.reinit(locally_relevant_dofs);

        double initial_residual = 1;
        double residual = 1;

        unsigned int newton_iteration = 0;
        for (; newton_iteration < max_newton_iterations; ++newton_iteration)
        {
          if (newton_iteration < 2)
            make_newton_constraints(constraints, 
                                    newton_constraints,
                                    displacements,
                                    (newton_iteration == 0));

          // Assemble and solve the linear system.
          assemble_linear_system(dof_handler, 
                                 newton_constraints,
                                 displacements,
                                 system_matrix,
                                 system_rhs,
                                 true);

          residual = system_rhs.l2_norm();
          if (residual < 1e-50)
          {
            displacements = 0;
            this->get_pcout() << "done." << std::endl;
            return;
          }

          solve_linear_system(dof_handler, 
                              newton_constraints,
                              system_matrix,
                              system_rhs,
                              newton_update);

          if (newton_iteration == 0)
          {
            dist_displacements = displacements;
            dist_displacements.add(1, newton_update);
            displacements = dist_displacements;

            initial_residual = residual;
            if (initial_residual <= std::numeric_limits<double>::min())
              break;
          }
          else
          {
            // Execute a line search if the solution update doesn't decrease the norm
            // of the rhs enough.
            double test_residual = 0;
            double lambda        = 1;
            double alpha         = 1e-4;

            unsigned int line_search_iteration = 0;
            for (; line_search_iteration < max_line_search_iterations; ++line_search_iteration)
            {
              dist_displacements = displacements;
              dist_displacements.sadd(1, lambda, newton_update);
              test_displacements = dist_displacements;

              assemble_linear_system(dof_handler, 
                                     newton_constraints,
                                     test_displacements,
                                     system_matrix,
                                     system_rhs,
                                     false);

              test_residual = system_rhs.l2_norm();
              if (test_residual < (1.0 - alpha * lambda) * residual)
              {
                residual = test_residual;
                break;
              }
              else
              {
                lambda *= (2.0 / 3.0);
              }
            }
          }

          displacements = dist_displacements;

          if (residual / initial_residual < newton_solver_tolerance)
            break;
        }

        this->get_pcout() << "done with " << newton_iteration << " iterations." << std::endl;
      }


      template <int dim>
      void Hyperelasticity<dim>::declare_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Mesh deformation");
        {
          prm.enter_subsection("Mesh smoothing");
          {
            prm.enter_subsection("Hyperelasticity");
            {
              prm.declare_entry("Poisson's ratio", "0.4",
                                Patterns::Double(-1.0, 0.5),
                                "");
              prm.declare_entry("Linear solver tolerance", "1e-7",
                                Patterns::Double(0),
                                "");
              prm.declare_entry("Newton solver tolerance", "1e-5",
                                Patterns::Double(0),
                                "");
              prm.declare_entry("Max Newton iterations", "100",
                                Patterns::Integer(0),
                                "");
              prm.declare_entry("Max line search iterations", "5",
                                Patterns::Integer(0),
                                "");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      Hyperelasticity<dim>::parse_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Mesh deformation");
        {
          prm.enter_subsection("Mesh smoothing");
          {
            prm.enter_subsection("Hyperelasticity");
            {
              double nu = prm.get_double("Poisson's ratio");
              // The Lame constant (mu is supposed to be 1).
              lambda = 2.0 * nu / (1.0 - 2.0 * nu);

              linear_solver_tolerance = prm.get_double("Linear solver tolerance");
              newton_solver_tolerance = prm.get_double("Newton solver tolerance");
              max_newton_iterations = prm.get_integer("Max Newton iterations");
              max_line_search_iterations = prm.get_integer("Max line search iterations");
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
    namespace MeshSmoothing
    {
      ELASPECT_REGISTER_MESH_SMOOTHING_MODEL(Hyperelasticity,
                                         "hyperelasticity",
                                         "")
    }
  }
}
