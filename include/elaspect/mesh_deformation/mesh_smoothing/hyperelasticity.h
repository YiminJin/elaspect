#ifndef _elaspect_mesh_deformation_mesh_smoothing_hyperelasticity_h
#define _elaspect_mesh_deformation_mesh_smoothing_hyperelasticity_h

#include <elaspect/mesh_deformation/mesh_smoothing/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace MeshDeformation
  {
    namespace MeshSmoothing
    {
      template <int dim>
      class Hyperelasticity : public Interface<dim>, public SimulatorAccess<dim>
      {
        public:
          void
          compute_mesh_displacements(const DoFHandler<dim>           &dof_handler,
                                     const AffineConstraints<double> &constraints,
                                     TrilinosWrappers::MPI::Vector   &displacements) const override;

          static
          void
          declare_parameters(ParameterHandler &prm);

          void
          parse_parameters(ParameterHandler &prm) override;

        private:
          void 
          make_newton_constraints(const AffineConstraints<double>     &constraints,
                                  AffineConstraints<double>           &newton_constraints,
                                  const TrilinosWrappers::MPI::Vector &displacements,
                                  const bool                           use_inhomogeneity) const;

          void 
          assemble_linear_system(const DoFHandler<dim>               &dof_handler,
                                 const AffineConstraints<double>     &newton_constraints,
                                 const TrilinosWrappers::MPI::Vector &displacements,
                                 TrilinosWrappers::SparseMatrix      &system_matrix,
                                 TrilinosWrappers::MPI::Vector       &system_rhs,
                                 const bool                           reassemble_matrix) const;

          void
          solve_linear_system(const DoFHandler<dim>                &dof_handler,
                              const AffineConstraints<double>      &newton_constraints,
                              const TrilinosWrappers::SparseMatrix &tangent_matrix,
                              const TrilinosWrappers::MPI::Vector  &system_rhs,
                              TrilinosWrappers::MPI::Vector        &newton_update) const;

          double lambda;

          double linear_solver_tolerance;

          double newton_solver_tolerance;

          unsigned int max_newton_iterations;

          unsigned int max_line_search_iterations;
      };
    }
  }
}

#endif
