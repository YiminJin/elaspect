#ifndef _elaspect_parameters_h
#define _elaspect_parameters_h

#include <elaspect/material_model/average.h>

#include <deal.II/base/parameter_handler.h>

namespace elaspect
{
  using namespace dealii;

  // forward declaration:
  namespace GeometryModel
  {
    template <int dim>
    class Interface;
  }

  namespace ConstitutiveRelation
  {
    enum Component
    {
      elasticity        = 0,
      viscosity         = 1,
      plasticity        = 2,
      thermal_expansion = 4,
      pore_fluid        = 8
    };

    inline Component operator | (const Component c1,
                                 const Component c2)
    {
      return static_cast<Component>(static_cast<unsigned int>(c1) |
                                    static_cast<unsigned int>(c2));
    }

    inline Component &operator |= (Component &c1,
                                   const Component c2)
    {
      c1 = c1 | c2;
      return c1;
    }

    inline Component operator & (const Component c1,
                                 const Component c2)
    {
      return static_cast<Component>(static_cast<unsigned int>(c1) &
                                    static_cast<unsigned int>(c2));
    }

    inline Component &operator &= (Component &c1,
                                   const Component c2)
    {
      c1 = c1 & c2;
      return c1;
    }
  }

  template <int dim>
  struct Parameters
  {
    struct LinearSolver
    {
      enum Kind
      {
        CG,
        BiCGStab,
        GMRES,
        MUMPS
      };
    };

    struct NullspaceRemoval
    {
      enum Kind
      {
        none = 0,
        net_translation_x = 0x1,
        net_translation_y = 0x2,
        net_translation_z = 0x4,
        net_translation   = 0x1+0x2+0x4,
        linear_momentum_x = 0x8,
        linear_momentum_y = 0x10,
        linear_momentum_z = 0x20,
        linear_momentum   = 0x8+0x10+0x20,
        net_rotation      = 0x40,
        angular_momentum  = 0x80
      };
    };

    struct HeatTransport
    {
      enum Kind
      {
        none,
        prescribed,
        convection_diffusion
      };
    };

    Parameters (ParameterHandler &prm,
                const MPI_Comm mpi_communicator);

    static
    void declare_parameters (ParameterHandler &prm);

    void parse_parameters (ParameterHandler &prm,
                           const MPI_Comm mpi_communicator);

    void parse_geometry_dependent_parameters (ParameterHandler &prm,
                                              const GeometryModel::Interface<dim> &geometry_model);

    typename ConstitutiveRelation::Component constitutive_relation;

    typename NullspaceRemoval::Kind nullspace_removal;
    typename HeatTransport::Kind    heat_transport;

    MaterialModel::MaterialAveraging::AveragingOperation material_averaging;

    std::string                 output_directory;
    bool                        convert_to_years;
    double                      start_time;
    unsigned int                timing_output_frequency;
    bool                        use_ALE_method;
    double                      CFL_number;
    bool                        run_postprocessors_on_initial_refinement;
    bool                        run_postprocessors_on_nonlinear_iterations;

    unsigned int                n_compositional_fields;
    std::vector<std::string>    names_of_compositional_fields;

    unsigned int                max_nonlinear_iterations;
    double                      nonlinear_tolerance;

    typename LinearSolver::Kind mechanical_system_linear_solver;
    double                      mechanical_system_solver_tolerance;
    unsigned int                mechanical_system_max_linear_iterations;
    unsigned int                mechanical_system_gmres_restart_length;
    double                      mechanical_system_amg_aggregation_threshold;
    unsigned int                max_line_search_steps;
    double                      line_search_beta;
    double                      line_search_delta;
    bool                        enforce_convergence_for_mechanical_system;

    double                      thermo_system_solver_tolerance;
    unsigned int                thermo_system_max_linear_iterations;

    unsigned int                initial_global_refinement;
    unsigned int                initial_adaptive_refinement;
    double                      refinement_fraction;
    double                      coarsening_fraction;
    unsigned int                min_grid_level;
    bool                        adapt_by_fraction_of_cells;
    unsigned int                adaptive_refinement_interval;
    bool                        skip_solvers_on_initial_refinement;

    unsigned int                displacement_degree;
    unsigned int                fluid_pressure_degree;
    unsigned int                temperature_degree;
    unsigned int                composition_degree;
    unsigned int                n_gaussian_points;

    bool                        apply_BP_limiter_to_compositional_fields;
    std::vector<double>         composition_max_preset;
    std::vector<double>         composition_min_preset;
    bool                        apply_WENO_limiter_to_physical_fields;
    double                      KXRCF_indicator_threshold;

    std::map<types::boundary_id, std::pair<std::string,std::string> > prescribed_traction_boundary_indicators;

    std::set<types::boundary_id> fixed_heat_flux_boundary_indicators;
  };
}

#endif
