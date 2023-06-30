#include <elaspect/simulator.h>

namespace elaspect
{
  template <int dim>
  SimulatorAccess<dim>::SimulatorAccess ()
    : simulator (nullptr)
  {}


  template <int dim>
  SimulatorAccess<dim>::SimulatorAccess (const Simulator<dim> &simulator_object)
    : simulator (&simulator_object)
  {}


  template <int dim>
  SimulatorAccess<dim>::~SimulatorAccess ()
  {}


  template <int dim>
  void
  SimulatorAccess<dim>::initialize_simulator (const Simulator<dim> &simulator_object)
  {
    simulator = &simulator_object;
  }


  template <int dim>
  const Simulator<dim> &
  SimulatorAccess<dim>::get_simulator () const
  {
    return *simulator;
  }


  template <int dim>
  const Introspection<dim> &
  SimulatorAccess<dim>::introspection () const
  {
    return simulator->introspection;
  }


  template <int dim>
  const Parameters<dim> &
  SimulatorAccess<dim>::get_parameters() const
  {
    return simulator->parameters;
  }


  template <int dim>
  SimulatorSignals<dim> &
  SimulatorAccess<dim>::get_signals() const
  {
    // Our reference to the Simulator is const, but we need to
    // be able to connect to the signals so a cast is required.
    return const_cast<SimulatorSignals<dim>&>(simulator->signals);
  }


  template <int dim>
  MPI_Comm
  SimulatorAccess<dim>::get_mpi_communicator () const
  {
    return simulator->mpi_communicator;
  }


  template <int dim>
  const ConditionalOStream &
  SimulatorAccess<dim>::get_pcout () const
  {
    return simulator->pcout;
  }


  template <int dim>
  TimerOutput &
  SimulatorAccess<dim>::get_computing_timer () const
  {
    return simulator->computing_timer;
  }


  template <int dim>
  double SimulatorAccess<dim>::get_time () const
  {
    return simulator->time;
  }


  template <int dim>
  double SimulatorAccess<dim>::get_timestep () const
  {
    return simulator->time_step;
  }


  template <int dim>
  unsigned int
  SimulatorAccess<dim>::get_timestep_number () const
  {
    return simulator->timestep_number;
  }


  template <int dim>
  unsigned int 
  SimulatorAccess<dim>::get_nonlinear_iteration () const
  {
    return simulator->nonlinear_iteration;
  }


  template <int dim>
  bool
  SimulatorAccess<dim>::convert_output_to_years () const
  {
    return simulator->parameters.convert_to_years;
  }


  template <int dim>
  unsigned int
  SimulatorAccess<dim>::n_compositional_fields () const
  {
    return simulator->parameters.n_compositional_fields;
  }


  template <int dim>
  ConstitutiveRelation::Component
  SimulatorAccess<dim>::constitutive_relation () const
  {
    return simulator->parameters.constitutive_relation;
  }


  template <int dim>
  unsigned int
  SimulatorAccess<dim>::get_displacement_degree () const
  {
    return simulator->parameters.displacement_degree;
  }


  template <int dim>
  std::string
  SimulatorAccess<dim>::get_output_directory () const
  {
    return simulator->parameters.output_directory;
  }


  template <int dim>
  const parallel::distributed::Triangulation<dim> &
  SimulatorAccess<dim>::get_triangulation () const
  {
    return simulator->triangulation;
  }


  template <int dim>
  double
  SimulatorAccess<dim>::get_volume () const
  {
    return simulator->global_volume;
  }


  template <int dim>
  const Mapping<dim> &
  SimulatorAccess<dim>::get_mapping () const
  {
    return *(simulator->mapping);
  }


  template <int dim>
  const FiniteElement<dim> &
  SimulatorAccess<dim>::get_fe () const
  {
    Assert (simulator->dof_handler.n_dofs() > 0,
            ExcMessage("You are trying to access the FiniteElement before the DOFs have been "
                       "initialized. This may happen when accessing the Simulator from a plugin "
                       "that gets executed early in some cases (like material models) or from "
                       "an early point in the core code."));
    return simulator->dof_handler.get_fe();
  }


  template <int dim>
  const DoFHandler<dim> &
  SimulatorAccess<dim>::get_dof_handler () const
  {
    return simulator->dof_handler;
  }


  template <int dim>
  const Quadrature<dim> &
  SimulatorAccess<dim>::get_quadrature_formula () const
  {
    return simulator->quadrature_formula;
  }


  template <int dim>
  const QPDHandler<dim> &
  SimulatorAccess<dim>::get_qpd_handler() const
  {
    return simulator->qpd_handler;
  }


  template <int dim>
  const TrilinosWrappers::MPI::BlockVector &
  SimulatorAccess<dim>::get_solution () const
  {
    return simulator->solution;
  }


  template <int dim>
  const TrilinosWrappers::MPI::BlockVector &
  SimulatorAccess<dim>::get_old_solution () const
  {
    return simulator->old_solution;
  }


  template <int dim>
  const GeometryModel::Interface<dim> &
  SimulatorAccess<dim>::get_geometry_model () const
  {
    Assert (simulator->geometry_model.get() != nullptr,
            ExcMessage("You can not call this function if no such model is actually available."));
    return *simulator->geometry_model.get();
  }


  template <int dim>
  const InitialTopography::Interface<dim> &
  SimulatorAccess<dim>::get_initial_topography () const
  {
    Assert (simulator->initial_topography.get() != nullptr,
            ExcMessage("You can not call this function if no such model is actually available."));
    return *simulator->initial_topography.get();
  }


  template <int dim>
  const HeatingModel::Manager<dim> &
  SimulatorAccess<dim>::get_heating_model_manager() const
  {
    return simulator->heating_model_manager;
  }


  template <int dim>
  const GravityModel::Interface<dim> &
  SimulatorAccess<dim>::get_gravity_model () const
  {
    Assert (simulator->gravity_model.get() != nullptr,
            ExcMessage("You can not call this function if no such model is actually available."));
    return *simulator->gravity_model.get();
  }


  template <int dim>
  const MaterialHandler<dim> &
  SimulatorAccess<dim>::get_material_handler () const
  {
    return simulator->material_handler;
  }


  template <int dim>
  const InitialComposition::Manager<dim> &
  SimulatorAccess<dim>::get_initial_composition_manager () const
  {
    return simulator->initial_composition_manager;
  }


  template <int dim>
  const InitialTemperature::Manager<dim> &
  SimulatorAccess<dim>::get_initial_temperature_manager () const
  {
    return simulator->initial_temperature_manager;
  }


  template <int dim>
  const BoundaryComposition::Manager<dim> &
  SimulatorAccess<dim>::get_boundary_composition_manager () const
  {
    return simulator->boundary_composition_manager;
  }


  template <int dim>
  const BoundaryTemperature::Manager<dim> &
  SimulatorAccess<dim>::get_boundary_temperature_manager () const
  {
    return simulator->boundary_temperature_manager;
  }


  template <int dim>
  const BoundaryHeatFlux::Interface<dim> &
  SimulatorAccess<dim>::get_boundary_heat_flux () const
  { 
    Assert (simulator->boundary_heat_flux.get() != nullptr,
            ExcMessage("You can not call this function if no such model is actually available."));
    return *simulator->boundary_heat_flux.get();
  }


  template <int dim>
  const BoundaryVelocity::Manager<dim> &
  SimulatorAccess<dim>::get_boundary_velocity_manager () const
  {
    return simulator->boundary_velocity_manager;
  }


  template <int dim>
  const std::map<types::boundary_id,std::unique_ptr<BoundaryTraction::Interface<dim>>> &
  SimulatorAccess<dim>::get_boundary_traction () const
  {
    return simulator->boundary_traction;
  }


  template <int dim>
  const std::set<types::boundary_id> &
  SimulatorAccess<dim>::get_fixed_composition_boundary_indicators () const
  {
    return get_boundary_composition_manager().get_fixed_composition_boundary_indicators();
  }


  template <int dim>
  const std::set<types::boundary_id> &
  SimulatorAccess<dim>::get_fixed_temperature_boundary_indicators () const
  {
    return get_boundary_temperature_manager().get_fixed_temperature_boundary_indicators();
  }


  template <int dim>
  const std::set<types::boundary_id> &
  SimulatorAccess<dim>::get_fixed_heat_flux_boundary_indicators () const
  {
    return simulator->parameters.fixed_heat_flux_boundary_indicators;
  }


  template <int dim>
  const Postprocess::Manager<dim> &
  SimulatorAccess<dim>::get_postprocess_manager() const
  {
    return simulator->postprocess_manager;
  }


  template <int dim>
  const MeshDeformationHandler<dim> &
  SimulatorAccess<dim>::get_mesh_deformation_handler() const
  {
    return simulator->mesh_deformation_handler;
  }


  template <int dim>
  bool
  SimulatorAccess<dim>::rebuild_sparsity_and_matrices() const
  {
    return simulator->rebuild_sparsity_and_matrices;
  }


  template <int dim>
  bool
  SimulatorAccess<dim>::reassemble_tangent_matrix() const
  {
    return simulator->reassemble_tangent_matrix;
  }
}


// explicit instantiations
namespace elaspect
{
#define INSTANTIATE(dim) \
  template class SimulatorAccess<dim>;

  ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
