#ifndef _elaspect_simulator_h
#define _elaspect_simulator_h

#include <elaspect/global.h>
#include <elaspect/quadrature_point_data.h>
#include <elaspect/simulator_access.h>
#include <elaspect/simulator_signals.h>
#include <elaspect/parameters.h>
#include <elaspect/introspection.h>
#include <elaspect/postprocess/interface.h>
#include <elaspect/mesh_refinement/interface.h>
#include <elaspect/geometry_model/interface.h>
#include <elaspect/initial_topography/interface.h>
#include <elaspect/mesh_deformation/handler.h>
#include <elaspect/material_model/handler.h>
#include <elaspect/heating_model/interface.h>
#include <elaspect/gravity_model/interface.h>
#include <elaspect/time_stepping/interface.h>
#include <elaspect/initial_composition/interface.h>
#include <elaspect/initial_temperature/interface.h>
#include <elaspect/boundary_composition/interface.h>
#include <elaspect/boundary_temperature/interface.h>
#include <elaspect/boundary_velocity/interface.h>
#include <elaspect/boundary_traction/interface.h>
#include <elaspect/boundary_heat_flux/interface.h>

#include <deal.II/base/timer.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1_eulerian.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>

#include <memory>
#include <thread>

namespace elaspect
{
  using namespace dealii;

  namespace internal
  {
    namespace Assembly
    {
      namespace Scratch
      {
        template <int dim> struct MechanicalSystem;
        template <int dim> struct HydroSystem;
        template <int dim> struct ThermoSystem;
        template <int dim> struct QPDSystem;
      }

      namespace CopyData
      {
        template <int dim> struct MechanicalSystem;
        template <int dim> struct HydroSystem;
        template <int dim> struct ThermoSystem;
        template <int dim> struct QPDSystem;
      }
    }
  }

  namespace Assemblers
  {
    template <int dim> class Interface;
    template <int dim> class Manager;
  }

  template <int dim>
  struct RotationProperties
  {
    RotationProperties()
      :
      scalar_moment_of_inertia(numbers::signaling_nan<double>()),
      scalar_angular_momentum(numbers::signaling_nan<double>()),
      scalar_rotation(numbers::signaling_nan<double>()),
      tensor_moment_of_inertia(numbers::signaling_nan<SymmetricTensor<2,dim>>()),
      tensor_angular_momentum(numbers::signaling_nan<Tensor<1,dim>>()),
      tensor_rotation(numbers::signaling_nan<Tensor<1,dim>>())
    {}

    /**
     * Scalar properties for the two-dimensional case
     * with a fixed rotation axis (z).
     */
    double scalar_moment_of_inertia;
    double scalar_angular_momentum;
    double scalar_rotation;

    /**
     * Tensor properties for the three-dimensional case.
     */
    SymmetricTensor<2,dim> tensor_moment_of_inertia;
    Tensor<1,dim> tensor_angular_momentum;
    Tensor<1,dim> tensor_rotation;
  };


  template <int dim>
  class Simulator
  {
    public:
      using LinearSolver     = typename Parameters<dim>::LinearSolver;
      using NullspaceRemoval = typename Parameters<dim>::NullspaceRemoval;

      Simulator (const MPI_Comm mpi_communicator,
                 ParameterHandler &prm);

      ~Simulator ();

      static
      void declare_parameters (ParameterHandler &prm);

      void run ();

      struct QPDField
      {
        enum FieldType
        {
          compositional_field,
          physical_field
        };

        FieldType field_type;

        const unsigned int component_index;

        const unsigned int block_index;

        const unsigned int base_index;

        const unsigned int qpd_indicator;

        const unsigned int compositional_variable;

        const FEValuesExtractors::Scalar scalar_extractor;

        QPDField (const Introspection<dim> &introspection,
                        const std::string        &name,
                        const unsigned int        component);

        bool is_compositional_field() const;
      };

    private:
      void setup_dofs ();

      void setup_introspection ();

      void setup_system_matrix ();

      void compute_current_constraints ();

      void 
      compute_current_displacement_boundary_constraints (AffineConstraints<double> &constraints);

      void initialize_temperature_field ();

      void initialize_fluid_pressure_field ();

      void initialize_compositional_fields ();

      void update_quadrature_point_data ();

      void setup_nullspace_constraints(AffineConstraints<double> &constraints);

      void setup_nullspace_constraints_for_sphere(AffineConstraints<double> &constraints);

      void remove_nullspace(TrilinosWrappers::MPI::BlockVector &du,
                            TrilinosWrappers::MPI::BlockVector &dist_du);

      RotationProperties<dim>
      compute_net_angular_momentum(const bool use_constant_density,
                                   const TrilinosWrappers::MPI::BlockVector &solution) const;

      void remove_net_linear_momentum( const bool use_constant_density,
                                       TrilinosWrappers::MPI::BlockVector &du,
                                       TrilinosWrappers::MPI::BlockVector &dist_du);

      void remove_net_angular_momentum( const bool use_constant_density,
                                        TrilinosWrappers::MPI::BlockVector &du,
                                        TrilinosWrappers::MPI::BlockVector &dist_du);

      void compute_initial_displacement_boundary_constraints (AffineConstraints<double> &constraints);

      void start_timestep ();

      bool assemble_and_solve_mechanical_system ();

      void assemble_and_solve_thermo_system ();

      void assemble_and_solve_hydro_system ();

      void assemble_and_solve_qpd_system ();

      void apply_return_mapping ();

      void postprocess ();

      void output_statistics ();

      void refine_mesh (const unsigned int max_grid_level);

      void set_assemblers ();

      void set_mechanical_system_assemblers ();

      void assemble_mechanical_system ();

      void 
      local_assemble_mechanical_system (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                   internal::Assembly::Scratch::MechanicalSystem<dim> &scratch,
                                   internal::Assembly::CopyData::MechanicalSystem<dim> &data);

      void
      copy_local_to_global_mechanical_system (const internal::Assembly::CopyData::MechanicalSystem<dim> &data);

      void solve_mechanical_system (TrilinosWrappers::MPI::BlockVector &dist_du);

      void assemble_thermo_system ();

      void
      local_assemble_thermo_system (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                     internal::Assembly::Scratch::ThermoSystem<dim> &scratch,
                                     internal::Assembly::CopyData::ThermoSystem<dim> &data);

      void
      copy_local_to_global_thermo_system (const internal::Assembly::CopyData::ThermoSystem<dim> &data);

      void solve_thermo_system ();

      void assemble_qpd_system (const std::vector<QPDField> &qpd_fields);

      void
      local_assemble_qpd_system (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                 internal::Assembly::Scratch::QPDSystem<dim> &scratch,
                                 internal::Assembly::CopyData::QPDSystem<dim> &data);

      void 
      copy_local_to_global_qpd_system (const internal::Assembly::CopyData::QPDSystem<dim> &data);

      void solve_qpd_system (const std::vector<QPDField> &qpd_fields);

      bool maybe_do_initial_refinement (const unsigned int max_refinement_level);

      void maybe_refine_mesh (const unsigned int max_refinement_level);

      void maybe_write_timing_output () const;

      void 
      interpolate_onto_displacement_system(const TensorFunction<1,dim> &func,
                                           TrilinosWrappers::MPI::Vector &vec);

      void interpolate_material_output_into_temperature_field ();

      void advect_quadrature_point_data();

      void refresh_quadrature_point_data();

      void apply_BP_limiter(const QPDField &qpd_field);

      template <typename T>
      void compute_KXRCF_indicators(Vector<T> &KXRCF_indicators,
                                    const QPDField &qpd_field) const;

      void apply_WENO_limiter(const QPDField &qpd_field);

      std::unique_ptr<Assemblers::Manager<dim>> assemblers;

      Parameters<dim>                           parameters;

      SimulatorSignals<dim>                     signals;

      Introspection<dim>                        introspection;

      MPI_Comm                                  mpi_communicator;

      std::ofstream log_file_stream;

      using TeeDevice = boost::iostreams::tee_device<std::ostream, std::ofstream>;
      using TeeStream = boost::iostreams::stream< TeeDevice >;

      TeeDevice iostream_tee_device;
      TeeStream iostream_tee_stream;

      ConditionalOStream    pcout;

      TableHandler          statistics;
      std::thread           output_statistics_thread;
      std::size_t           statistics_last_write_size;
      std::size_t           statistics_last_hash;

      mutable TimerOutput   computing_timer;

      const std::unique_ptr<GeometryModel::Interface<dim>>      geometry_model;
      const std::unique_ptr<InitialTopography::Interface<dim>>  initial_topography;
      const std::unique_ptr<GravityModel::Interface<dim>>       gravity_model;
      const std::unique_ptr<BoundaryHeatFlux::Interface<dim>>   boundary_heat_flux;
      std::map<types::boundary_id,std::unique_ptr<BoundaryTraction::Interface<dim>>> boundary_traction;

      MaterialHandler<dim>                                  material_handler;

      HeatingModel::Manager<dim>                            heating_model_manager;
      TimeStepping::Manager<dim>                            time_stepping_manager;
      InitialTemperature::Manager<dim>                      initial_temperature_manager;
      InitialComposition::Manager<dim>                      initial_composition_manager;
      BoundaryTemperature::Manager<dim>                     boundary_temperature_manager;
      BoundaryComposition::Manager<dim>                     boundary_composition_manager;
      BoundaryVelocity::Manager<dim>                        boundary_velocity_manager;
      MeshRefinement::Manager<dim>                          mesh_refinement_manager;
      Postprocess::Manager<dim>                             postprocess_manager;

      double                time;
      double                time_step;
      unsigned int          timestep_number;
      unsigned int          pre_refinement_step;
      unsigned int          nonlinear_iteration;

      parallel::distributed::Triangulation<dim>             triangulation;
      double                                                global_volume;

      MeshDeformationHandler<dim>                           mesh_deformation_handler;

      const FESystem<dim>                                   finite_element;

      DoFHandler<dim>                                       dof_handler;

      std::unique_ptr<MappingQ1Eulerian<dim, TrilinosWrappers::MPI::Vector>> mapping;

      QGauss<dim>                                           quadrature_formula;

      QPDHandler<dim>                                       qpd_handler;

      AffineConstraints<double>                             constraints;
      AffineConstraints<double>                             current_constraints;

      TrilinosWrappers::BlockSparseMatrix                   system_matrix;

      TrilinosWrappers::MPI::BlockVector                    solution;
      TrilinosWrappers::MPI::BlockVector                    old_solution;
      
      TrilinosWrappers::MPI::BlockVector                    system_rhs;

      bool rebuild_sparsity_and_matrices;
      bool reassemble_tangent_matrix;

      friend class SimulatorAccess<dim>;
  };
}

#endif
