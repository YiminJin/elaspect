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


#ifndef _elaspect_simulator_access_h
#define _elaspect_simulator_access_h

#include <elaspect/global.h>
#include <elaspect/parameters.h>
#include <elaspect/introspection.h>

#include <deal.II/base/timer.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>

namespace elaspect
{
  using namespace dealii;

  template <int dim> class Simulator;
  template <int dim> struct SimulatorSignals;
  
  namespace GeometryModel
  {
    template <int dim> class Interface;
  }

  namespace InitialTopography
  {
    template <int dim> class Interface;
  }

  namespace HeatingModel
  {
    template <int dim> class Manager;
  }

  namespace MaterialModel
  {
    template <int dim> class Interface;
    template <int dim> struct MaterialModelInputs;
    template <int dim> struct MaterialModelOutputs;
  }

  namespace GravityModel
  {
    template <int dim> class Interface;
  }

  namespace InitialComposition
  {
    template <int dim> class Manager;
  }

  namespace InitialTemperature
  {
    template <int dim> class Manager;
  }

  namespace BoundaryComposition
  {
    template <int dim> class Manager;
  }

  namespace BoundaryTemperature
  {
    template <int dim> class Manager;
  }

  namespace BoundaryHeatFlux
  {
    template <int dim> class Interface;
  }

  namespace BoundaryVelocity
  {
    template <int dim> class Manager;
  }

  namespace BoundaryTraction
  {
    template <int dim> class Interface;
  }

  namespace Postprocess
  {
    template <int dim> class Manager;
  }

  template <int dim> class MaterialHandler;
  template <int dim> class MeshDeformationHandler;

  template <int dim>
  class SimulatorAccess
  {
    public:
      SimulatorAccess ();

      SimulatorAccess (const Simulator<dim> &simulator_object);

      virtual
      ~SimulatorAccess ();

      virtual
      void
      initialize_simulator (const Simulator<dim> &simulator_object);

      const Introspection<dim> &
      introspection () const;

      const Simulator<dim> &
      get_simulator () const;

      const Parameters<dim> &
      get_parameters () const;

      SimulatorSignals<dim> &
      get_signals() const;

      MPI_Comm
      get_mpi_communicator () const;

      const ConditionalOStream &
      get_pcout () const;

      TimerOutput &
      get_computing_timer () const;

      double get_time () const;

      double get_timestep () const;

      unsigned int
      get_timestep_number () const;

      unsigned int
      get_nonlinear_iteration () const;

      bool
      convert_output_to_years () const;

      std::string
      get_output_directory () const;

      unsigned int
      n_compositional_fields () const;

      ConstitutiveRelation::Component
      constitutive_relation () const;

      unsigned int 
      get_displacement_degree () const;

      const parallel::distributed::Triangulation<dim> &
      get_triangulation () const;

      double
      get_volume () const;
      
      const Mapping<dim> &
      get_mapping () const;

      const FiniteElement<dim> &
      get_fe () const;

      const DoFHandler<dim> &
      get_dof_handler () const;

      const Quadrature<dim> &
      get_quadrature_formula () const;

      const QPDHandler<dim> &
      get_qpd_handler () const;

      const TrilinosWrappers::MPI::BlockVector &
      get_solution () const;

      const TrilinosWrappers::MPI::BlockVector &
      get_old_solution () const;

      const GeometryModel::Interface<dim> &
      get_geometry_model () const;

      const InitialTopography::Interface<dim> &
      get_initial_topography () const;

      const HeatingModel::Manager<dim> &
      get_heating_model_manager() const;

      const GravityModel::Interface<dim> &
      get_gravity_model () const;

      const MaterialHandler<dim> &
      get_material_handler () const;

      const InitialComposition::Manager<dim> &
      get_initial_composition_manager () const;

      const InitialTemperature::Manager<dim> &
      get_initial_temperature_manager () const;

      const BoundaryComposition::Manager<dim> &
      get_boundary_composition_manager () const;

      const BoundaryTemperature::Manager<dim> &
      get_boundary_temperature_manager () const;

      const BoundaryHeatFlux::Interface<dim> &
      get_boundary_heat_flux () const;

      const BoundaryVelocity::Manager<dim> &
      get_boundary_velocity_manager () const;

      const std::map<types::boundary_id, std::unique_ptr<BoundaryTraction::Interface<dim>>> &
      get_boundary_traction () const;

      const std::set<types::boundary_id> &
      get_fixed_composition_boundary_indicators () const;

      const std::set<types::boundary_id> &
      get_fixed_temperature_boundary_indicators () const;

      const std::set<types::boundary_id> &
      get_fixed_heat_flux_boundary_indicators () const;

      const Postprocess::Manager<dim> &
      get_postprocess_manager() const;

      const MeshDeformationHandler<dim> &
      get_mesh_deformation_handler() const;

      bool
      rebuild_sparsity_and_matrices() const;

      bool
      reassemble_tangent_matrix() const;

    private:
      const Simulator<dim> *simulator;
  };
}

#endif
