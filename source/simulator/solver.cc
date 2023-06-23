#include <elaspect/simulator.h>

#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_precondition.h>

namespace elaspect
{
  namespace
  {
    void iterative_solver_failed(const std::string &solver_name,
                                 const std::string &output_filename,
                                 const std::vector<SolverControl> &solver_controls,
                                 const std::exception &exc)
    {
      // output solver history
      std::ofstream f((output_filename).c_str());

      for (unsigned int i=0; i<solver_controls.size(); ++i)
        {
          if (i>0)
            f << "\n";

          // Only request the solver history if a history has actually been created
          for (unsigned int j=0; j<solver_controls[i].get_history_data().size(); ++j)
            f << j << " " << solver_controls[i].get_history_data()[j] << "\n";
        }

      f.close();

      AssertThrow (false,
                   ExcMessage ("The " + solver_name
                               + " did not converge. It reported the following error:\n\n"
                               +
                               exc.what()
                               + "\n The required residual for convergence is: " 
                               + std::to_string(solver_controls.front().tolerance())
                               + ".\n See " + output_filename
                               + " for convergence history."));
    }
  }


  template <int dim>
  void 
  Simulator<dim>::
  solve_mechanical_system (TrilinosWrappers::MPI::BlockVector &distributed_solution_vector)
  {
    TimerOutput::Scope timer (computing_timer, "Solve mechanical system");

    const unsigned int block_idx = introspection.block_indices.displacement;

    // check if RHS is zero
    const double rhs_norm = system_rhs.block(block_idx).l2_norm();
    if (rhs_norm <= std::numeric_limits<double>::min())
    {
      pcout << "   Skip solving mechanical system because RHS is zero." << std::endl;
      distributed_solution_vector.block(block_idx) = 0;
      return;
    }

    if (parameters.mechanical_system_linear_solver == LinearSolver::MUMPS)
    {
      pcout << "   Solving mechanical system with MUMPS solver... " << std::flush;

      SolverControl solver_control;
      TrilinosWrappers::SolverDirect::AdditionalData solver_settings(
        false, "Amesos_Mumps");
      TrilinosWrappers::SolverDirect solver(solver_control, solver_settings);

      try
      {
        solver.solve(system_matrix.block(block_idx, block_idx),
                     distributed_solution_vector.block(block_idx),
                     system_rhs.block(block_idx));
      }
      catch (const std::exception &exc)
      {
        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
        {
          AssertThrow(false,
                      ExcMessage("The MUMPS solver failed to solve the linear system. "
                                 "The error code is:\n" + std::string(exc.what())));
        }
        else
          throw QuietException();
      }

      pcout << "done." << std::endl;
    }
    else
    {
      const double tolerance = std::max(1e-50, parameters.mechanical_system_solver_tolerance * rhs_norm);
      SolverControl solver_control(parameters.mechanical_system_max_linear_iterations, tolerance);
      solver_control.enable_history_data();

      std::unique_ptr<TrilinosWrappers::SolverBase> solver;
      std::string solver_name;

      switch (parameters.mechanical_system_linear_solver)
      {
        case LinearSolver::CG:
        {
          solver = std::make_unique<TrilinosWrappers::SolverCG>(solver_control);
          solver_name = "CG";
          break;
        }
        case LinearSolver::BiCGStab:
        {
          solver = std::make_unique<TrilinosWrappers::SolverBicgstab>(solver_control);
          solver_name = "BiCGStab";
          break;
        }
        case LinearSolver::GMRES:
        {
          TrilinosWrappers::SolverGMRES::AdditionalData solver_settings(
            false, parameters.mechanical_system_gmres_restart_length);
          solver = std::make_unique<TrilinosWrappers::SolverGMRES>(solver_control,
                                                                   solver_settings);
          solver_name = "GMRES";
          break;
        }
        default:
          Assert(false, ExcInternalError());
      }

      pcout << "   Solving mechanical system with " << solver_name << " solver... " << std::flush;

      TrilinosWrappers::PreconditionAMG preconditioner;
      TrilinosWrappers::PreconditionAMG::AdditionalData Amg_data;

      std::vector<std::vector<bool> > constant_modes;
      DoFTools::extract_constant_modes (dof_handler,
                                        introspection.component_masks.displacement,
                                        constant_modes);
      Amg_data.constant_modes = constant_modes;
      Amg_data.elliptic = true;
      Amg_data.higher_order_elements = (parameters.displacement_degree > 1);
      Amg_data.n_cycles = 1;
      Amg_data.w_cycle  = false;
      Amg_data.smoother_sweeps = 2;
      Amg_data.aggregation_threshold = parameters.mechanical_system_amg_aggregation_threshold;

      preconditioner.initialize(system_matrix.block(block_idx,block_idx), Amg_data);

     try
      {
        solver->solve(system_matrix.block(block_idx, block_idx),
                      distributed_solution_vector.block(block_idx),
                      system_rhs.block(block_idx),
                      preconditioner);
      }
      catch (const std::exception &exc)
      {
        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
        {
          iterative_solver_failed(solver_name + " solver for mechanical system",
                                  parameters.output_directory+"solver_history.txt",
                                  std::vector<SolverControl> {solver_control},
                                  exc);
        }
        else
          throw QuietException();
      }

      pcout << solver_control.last_step() << " iterations." << std::endl;
    }

    current_constraints.distribute(distributed_solution_vector);

    // remove nullspace
    TrilinosWrappers::MPI::BlockVector 
    solution_vector(introspection.index_sets.system_partitioning, 
                    introspection.index_sets.system_relevant_partitioning, 
                    mpi_communicator);
    solution_vector.block(block_idx) = distributed_solution_vector.block(block_idx);
    remove_nullspace (solution_vector, distributed_solution_vector);
  }


  template <int dim>
  void Simulator<dim>::solve_thermo_system ()
  {
    TimerOutput::Scope timer (computing_timer, "Solve thermal system");

    const unsigned int block_idx = introspection.block_indices.temperature;

    const double tolerance = std::max(1e-50, parameters.thermo_system_solver_tolerance * system_rhs.block(block_idx).l2_norm());
    SolverControl solver_control (parameters.thermo_system_max_linear_iterations, tolerance);
    solver_control.enable_history_data();

    std::unique_ptr<TrilinosWrappers::SolverBase> solver;
    std::string solver_name;
    if (parameters.use_ALE_method)
    {
      solver = std::make_unique<TrilinosWrappers::SolverGMRES>(solver_control);
      solver_name = "GMRES";
    }
    else
    {
      solver = std::make_unique<TrilinosWrappers::SolverCG>(solver_control);
      solver_name = "CG";
    }

    TrilinosWrappers::PreconditionILU preconditioner;
    preconditioner.initialize(system_matrix.block(block_idx, block_idx));

    pcout << "   Solving thermal system with " << solver_name 
          << " solver... " << std::flush;

    TrilinosWrappers::MPI::BlockVector
    distributed_solution (introspection.index_sets.system_partitioning, mpi_communicator);
    distributed_solution.block(block_idx) = solution.block(block_idx);
    current_constraints.set_zero (distributed_solution);

    try
    {
      solver->solve(system_matrix.block(block_idx, block_idx),
                    distributed_solution.block(block_idx),
                    system_rhs.block(block_idx),
                    preconditioner);
    }
    catch (const std::exception &exc)
    {
      if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        iterative_solver_failed("thermal system linear solver",
                                parameters.output_directory+"solver_history.txt",
                                std::vector<SolverControl> {solver_control},
                                exc);
      }
      else
        throw QuietException();
    }

    current_constraints.distribute (distributed_solution);
    solution.block(block_idx) = distributed_solution.block(block_idx);

    // print number of iterations and also record it in the
    // statistics file
    pcout << solver_control.last_step()
          << " iterations." << std::endl;
  }


  template <int dim>
  void
  Simulator<dim>::solve_qpd_system(const std::vector<QPDField> &qpd_fields)
  {
    if (qpd_fields.size() == 0)
      return;

    // Check if all the fields are in the same category.
#ifdef DEBUG
    const bool is_compositional_field = qpd_fields[0].is_compositional_field();
    for (unsigned int i = 1; i < qpd_fields.size(); ++i)
      Assert(qpd_fields[i].is_compositional_field() == is_compositional_field,
             ExcInternalError());
#endif

    // all QPD fields in the same category use the same block (the block of the
    // first field in the category) of system matrix for solving.
    const unsigned int block0_idx = qpd_fields[0].block_index;

    // build the preconditioner.
    TrilinosWrappers::PreconditionILU preconditioner;
    preconditioner.initialize(system_matrix.block(block0_idx, block0_idx));

    TrilinosWrappers::MPI::BlockVector distributed_solution(
      introspection.index_sets.system_partitioning,
      mpi_communicator);

    for (unsigned int field = 0; field < qpd_fields.size(); ++field)
    {
      // check if RHS is zero
      const unsigned int block_idx = qpd_fields[field].block_index;
      const double tolerance = 1e-8 * system_rhs.block(block_idx).l2_norm();
      if (tolerance <= std::numeric_limits<double>::min())
      {
        distributed_solution.block(block_idx) = 0;
        continue;
      }

      SolverControl solver_control(1000, tolerance);

      std::unique_ptr<TrilinosWrappers::SolverBase> solver;
      std::string solver_name;
      if (parameters.use_ALE_method)
      {
        solver = std::make_unique<TrilinosWrappers::SolverGMRES>(solver_control);
        solver_name = "GMRES";
      }
      else
      {
        solver = std::make_unique<TrilinosWrappers::SolverCG>(solver_control);
        solver_name = "CG";
      }

      // solve the linear system
      try
      {
        solver->solve(system_matrix.block(block0_idx, block0_idx),
                      distributed_solution.block(block_idx),
                      system_rhs.block(block_idx),
                      preconditioner);
      }
      // if the solver fails, report the error from processor 0 with some additional
      // information about its location, and throw a quiet exception on all other
      // processors
      catch (const std::exception &exc)
      {
        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
        {
          iterative_solver_failed(solver_name + " solver for QPD system",
                                  parameters.output_directory + "solver_history.txt",
                                  std::vector<SolverControl>{solver_control},
                                  exc);
        }
        else
          throw QuietException();
      }

      current_constraints.distribute(distributed_solution);
      solution.block(block_idx) = distributed_solution.block(block_idx);

      // apply BP/WENO limiter if requested
      if (qpd_fields[field].is_compositional_field())
      {
        if (parameters.apply_BP_limiter_to_compositional_fields)
          apply_BP_limiter(qpd_fields[field]);
      }
      else
      {
        if (parameters.apply_WENO_limiter_to_physical_fields)
          apply_WENO_limiter(qpd_fields[field]);
      }
    }
  }
}


// explicit instantiation
namespace elaspect
{
  #define INSTANTIATE(dim) \
  template void Simulator<dim>::solve_mechanical_system(TrilinosWrappers::MPI::BlockVector &); \
  template void Simulator<dim>::solve_thermo_system(); \
  template void Simulator<dim>::solve_qpd_system(const std::vector<QPDField> &);

  ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
