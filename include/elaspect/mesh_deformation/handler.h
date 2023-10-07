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


#ifndef _elaspect_mesh_deformation_handler_h
#define _elaspect_mesh_deformation_handler_h

#include <elaspect/global.h>
#include <elaspect/simulator.h>
#include <elaspect/mesh_deformation/free_surface/interface.h>
#include <elaspect/mesh_deformation/mesh_smoothing/interface.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/lac/block_vector.h>

namespace elaspect
{
  using namespace dealii;

  template <int dim> class SimulatorAccess;

  template <int dim>
  class MeshDeformationHandler : public SimulatorAccess<dim>
  {
    public:
      MeshDeformationHandler ();

      ~MeshDeformationHandler() override;

      void update();

      void execute();

      void setup_dofs();

      static
      void declare_parameters (ParameterHandler &prm);

      void parse_parameters (ParameterHandler &prm);

      const std::set<types::boundary_id> &
      get_free_surface_boundary_indicators() const;

      const TrilinosWrappers::MPI::BlockVector &
      get_mesh_displacement_increments() const;

    private:
      void make_constraints();

      void interpolate_mesh_displacement_increments();

      void interpolate_material_deformation_onto_vertices();

      void set_initial_topography ();

      std::unique_ptr<MeshDeformation::FreeSurface::Interface<dim>> free_surface_model;

      std::unique_ptr<MeshDeformation::MeshSmoothing::Interface<dim>> mesh_smoothing_model;

      const FESystem<dim> mesh_deformation_fe;

      DoFHandler<dim> mesh_deformation_dof_handler;

      TrilinosWrappers::MPI::BlockVector mesh_displacement_increments;

      TrilinosWrappers::MPI::Vector mesh_displacements;

      TrilinosWrappers::MPI::Vector old_mesh_displacements;

      TrilinosWrappers::MPI::Vector initial_topography;

      IndexSet mesh_locally_owned;

      IndexSet mesh_locally_relevant;

      AffineConstraints<double> mesh_deformation_constraints;

      AffineConstraints<double> mesh_vertex_constraints;

      std::set<types::boundary_id> tangential_mesh_deformation_boundary_indicators;

      std::set<types::boundary_id> zero_mesh_deformation_boundary_indicators;

      std::set<types::boundary_id> free_surface_boundary_indicators;

      bool include_initial_topography;

      double free_surface_theta;

      friend class Simulator<dim>;
  };
}

#endif
