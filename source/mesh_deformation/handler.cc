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


#include <elaspect/mesh_deformation/handler.h>
#include <elaspect/simulator_access.h>
#include <elaspect/geometry_model/interface.h>
#include <elaspect/initial_topography/zero_topography.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1_eulerian.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/numerics/vector_tools.h>

namespace elaspect
{
  template <int dim>
  MeshDeformationHandler<dim>::MeshDeformationHandler ()
    :
    mesh_deformation_fe (FE_Q<dim>(1), dim),
    include_initial_topography(false)
  {}


  template <int dim>
  MeshDeformationHandler<dim>::~MeshDeformationHandler()
  {}


  template <int dim>
  void MeshDeformationHandler<dim>::update ()
  {
    // Let the mesh smoothing model update itself.
    // (The free surface model will be updated in function execute().)
    mesh_smoothing_model->update();
  }


  template <int dim>
  void MeshDeformationHandler<dim>::execute ()
  {
    TimerOutput::Scope timer (this->get_computing_timer(), "Mesh deformation");

    if (this->get_parameters().use_ALE_method)
    {
      old_mesh_displacements = mesh_displacements;

      // Make the constraints for the mesh smoothing model.
      make_constraints();

      // Execute mesh smoothing
      mesh_smoothing_model->compute_mesh_displacements(mesh_deformation_dof_handler,
                                                       mesh_deformation_constraints,
                                                       mesh_displacements);

      // Interpolate the mesh displacement increments onto the same 
      // finite element space as used in the simulator,  which is 
      // needed for the ALE corrections.
      interpolate_mesh_displacement_increments();
    }
    else
    {
      // If ALE method is not applied, then directly interpolate 
      // the displacement solution onto mesh vertices.
      interpolate_material_deformation_onto_vertices();
    }
  }


  template <int dim>
  void MeshDeformationHandler<dim>::make_constraints ()
  {
    mesh_deformation_constraints.clear();
    mesh_deformation_constraints.reinit(mesh_locally_relevant);

    // mesh_velocity_constraints can use the same hanging node
    // information that was used for mesh_vertex constraints.
    mesh_deformation_constraints.merge(mesh_vertex_constraints);

    // Zero out the displacement for the zero-velocity boundaries
    // if the boundary is not in the set of tangential mesh boundaries 
    // and not in the set of mesh deformation boundary indicators.
    for (const auto &boundary_id : zero_mesh_deformation_boundary_indicators)
    {
      VectorTools::interpolate_boundary_values (this->get_mapping(),
                                                mesh_deformation_dof_handler,
                                                boundary_id,
                                                Functions::ZeroFunction<dim>(dim),
                                                mesh_deformation_constraints);
    }

    // Make the tangential deformation boundary constraints
    VectorTools::compute_no_normal_flux_constraints (mesh_deformation_dof_handler,
                                                     /* first_vector_component= */0,
                                                     tangential_mesh_deformation_boundary_indicators,
                                                     mesh_deformation_constraints,
                                                     this->get_mapping());

    // Ask the free surface model to make constraints for 
    // free surface boudaries.
    // Notice that the free surface model should be updated here
    // because the update function may depend on the solution of
    // incremental displacement.
    // TODO
    free_surface_model->update();

    AffineConstraints<double> free_surface_constraints(mesh_vertex_constraints.get_local_lines());
    free_surface_model->make_boundary_constraints(mesh_displacements,
                                                  mesh_deformation_dof_handler,
                                                  free_surface_constraints,
                                                  free_surface_boundary_indicators);

    mesh_deformation_constraints.merge(free_surface_constraints,
                                       AffineConstraints<double>::left_object_wins);
    mesh_deformation_constraints.close();
  }


  template <int dim>
  void MeshDeformationHandler<dim>::interpolate_material_deformation_onto_vertices()
  {
    mesh_deformation_constraints.clear();
    mesh_deformation_constraints.reinit(mesh_locally_relevant);

    // mesh_velocity_constraints can use the same hanging node
    // information that was used for mesh_vertex constraints.
    mesh_deformation_constraints.merge(mesh_vertex_constraints);
    mesh_deformation_constraints.close();

    std::vector<types::global_dof_index> vertex_dof_indices(dim), vertex_md_dof_indices(dim);
    std::vector<double> vertex_displacements(dim), incremental_displacements(dim);

    std::vector<bool> vertex_visited (this->get_triangulation().n_vertices(), false);
    TrilinosWrappers::MPI::Vector distributed_mesh_displacements(
      mesh_locally_owned, this->get_mpi_communicator());

    typename DoFHandler<dim>::active_cell_iterator
    cell = this->get_dof_handler().begin_active(),
    md_cell = mesh_deformation_dof_handler.begin_active();
    for (; cell != this->get_dof_handler().end(); ++cell, ++md_cell)
      if (cell->is_locally_owned())
      {
        for (const unsigned int vertex : cell->vertex_indices())
          if (!vertex_visited[cell->vertex_index(vertex)])
          {
            vertex_visited[cell->vertex_index(vertex)] = true;

            for (unsigned int d = 0; d < dim; ++d)
            {
              vertex_dof_indices[d]    = cell->vertex_dof_index (vertex, d);
              vertex_md_dof_indices[d] = md_cell->vertex_dof_index (vertex, d);
            }

            // skip the hanging nodes
            if (mesh_deformation_constraints.is_constrained(vertex_md_dof_indices[0]))
              continue;

            mesh_displacements.extract_subvector_to(vertex_md_dof_indices, vertex_displacements);
            this->get_solution().extract_subvector_to(vertex_dof_indices, incremental_displacements);

            for (unsigned int d = 0; d < dim; ++d)
              vertex_displacements[d] += incremental_displacements[d];

            distributed_mesh_displacements.set(vertex_md_dof_indices, vertex_displacements);
          }
      }

    distributed_mesh_displacements.compress(VectorOperation::insert);
    mesh_deformation_constraints.distribute(distributed_mesh_displacements);
    mesh_displacements = distributed_mesh_displacements;
  }


  template <int dim>
  void MeshDeformationHandler<dim>::setup_dofs ()
  {
    // these live in the same FE as the global system:
    mesh_displacement_increments.reinit(this->introspection().index_sets.system_partitioning,
                                        this->introspection().index_sets.system_relevant_partitioning,
                                        this->get_mpi_communicator());

    mesh_deformation_dof_handler.distribute_dofs (mesh_deformation_fe);

    // Renumber the DoFs hierarchical so that we get the
    // same numbering if we resume the computation. This
    // is because the numbering depends on the order the
    // cells are created.
    DoFRenumbering::hierarchical (mesh_deformation_dof_handler);

    this->get_pcout() << "Number of mesh deformation degrees of freedom: "
                      << mesh_deformation_dof_handler.n_dofs()
                      << std::endl;

    mesh_locally_owned = mesh_deformation_dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs (mesh_deformation_dof_handler,
                                             mesh_locally_relevant);

    // initialize the mesh displacement and free surface mesh velocity vectors
    mesh_displacements.reinit(mesh_locally_owned, mesh_locally_relevant, this->get_mpi_communicator());
    old_mesh_displacements.reinit(mesh_locally_owned, mesh_locally_relevant, this->get_mpi_communicator());
    initial_topography.reinit(mesh_locally_owned, mesh_locally_relevant, this->get_mpi_communicator());

    // if we are just starting, we need to set the initial topography.
    if (this->get_timestep_number() == 0)
      set_initial_topography();

    // We would like to make sure that the mesh stays conforming upon
    // redistribution, so we construct mesh_vertex_constraints, which
    // keeps track of hanging node constraints.
    // Note: this would be a more natural fit in make_constraints(),
    // but we would like to be able to apply vertex constraints directly
    // after setup_dofs(), as is done, for instance, during mesh
    // refinement.
    mesh_vertex_constraints.clear();
    mesh_vertex_constraints.reinit(mesh_locally_relevant);

    DoFTools::make_hanging_node_constraints(mesh_deformation_dof_handler, 
                                            mesh_vertex_constraints);

    mesh_vertex_constraints.close();

  }


  template <int dim>
  void MeshDeformationHandler<dim>::set_initial_topography ()
  {
    TrilinosWrappers::MPI::Vector distributed_initial_topography;
    distributed_initial_topography.reinit (mesh_locally_owned, this->get_mpi_communicator());

    if (!include_initial_topography)
      distributed_initial_topography = 0;
    else
    {
      const std::vector<Point<dim> > support_points 
        = mesh_deformation_fe.base_element(0).get_unit_support_points();

      const Quadrature<dim> quad(support_points);
      FEValues<dim> fe_values (this->get_mapping(), mesh_deformation_fe, quad, update_quadrature_points);

      const unsigned int n_q_points = fe_values.n_quadrature_points,
                         dofs_per_cell = fe_values.dofs_per_cell;

      std::vector<types::global_dof_index> cell_dof_indices (dofs_per_cell);

      for (const auto &cell : mesh_deformation_dof_handler.active_cell_iterators())
        if (cell->is_locally_owned())
        {
          cell->get_dof_indices (cell_dof_indices);

          fe_values.reinit (cell);
          for (unsigned int q = 0; q < n_q_points; ++q)
          {
            Point<dim-1> surface_point;
            unsigned int support_point_index = numbers::invalid_unsigned_int;
            std::array<double, dim> natural_coord = 
              this->get_geometry_model().cartesian_to_natural_coordinates(fe_values.quadrature_point(q));
            if (this->get_geometry_model().natural_coordinate_system()
                == Utilities::Coordinates::cartesian)
            {
              for (unsigned int d = 0; d < dim-1; ++d)
                surface_point[d] = natural_coord[d];
              support_point_index = mesh_deformation_fe.component_to_system_index (dim-1, q);
            }
            else if (this->get_geometry_model().natural_coordinate_system()
                     == Utilities::Coordinates::spherical)
            {
              for (unsigned int d = 1; d < dim; ++d)
                surface_point[d-1] = natural_coord[d];
              support_point_index = mesh_deformation_fe.component_to_system_index (0, q);
            }
            else
              AssertThrow (false, ExcNotImplemented());

            // Get the topography at this point.
            distributed_initial_topography[cell_dof_indices[support_point_index]]
              = this->get_initial_topography().value(surface_point);
          }
        }
    }

    distributed_initial_topography.compress (VectorOperation::insert);
    initial_topography = distributed_initial_topography;
  }


  template <int dim>
  void MeshDeformationHandler<dim>::interpolate_mesh_displacement_increments()
  {
    const Introspection<dim> &introspection = this->introspection();
    const FiniteElement<dim> &fe = this->get_fe();
    TrilinosWrappers::MPI::BlockVector distributed_vector(introspection.index_sets.system_partitioning,
                                                          this->get_mpi_communicator());

    const std::vector<Point<dim>> support_points =
      fe.base_element(introspection.base_elements.displacement).get_unit_support_points();

    const Quadrature<dim> quadrature(support_points);
    const UpdateFlags update_flags(update_values | update_JxW_values);
    FEValues<dim> md_fe_values(this->get_mapping(), mesh_deformation_fe, quadrature, update_flags);
    FEValues<dim> fe_values(this->get_mapping(), fe, quadrature, update_flags);

    const FEValuesExtractors::Vector extractor(0);

    const unsigned int n_q_points = fe_values.n_quadrature_points;
    const unsigned int dofs_per_cell = fe_values.dofs_per_cell;

    std::vector<types::global_dof_index> cell_dof_indices(dofs_per_cell);
    std::vector<Tensor<1,dim>> mesh_displacement_values(n_q_points);
    std::vector<Tensor<1,dim>> old_mesh_displacement_values(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator
    cell = this->get_dof_handler().begin_active(), 
    endc = this->get_dof_handler().end();
    typename DoFHandler<dim>::active_cell_iterator
    md_cell = mesh_deformation_dof_handler.begin_active();
    for (; cell != endc; ++cell, ++md_cell)
      if (cell->is_locally_owned())
      {
        cell->get_dof_indices(cell_dof_indices);

        fe_values.reinit(cell);
        md_fe_values.reinit(md_cell);
        md_fe_values[extractor].get_function_values(mesh_displacements, mesh_displacement_values);
        md_fe_values[extractor].get_function_values(old_mesh_displacements, old_mesh_displacement_values);

        for (unsigned int q = 0; q < n_q_points; ++q)
          for (unsigned int d = 0; d < dim; ++d)
          {
            const unsigned int support_point_index =
              fe.component_to_system_index(introspection.component_indices.displacement[d], q);
            distributed_vector[cell_dof_indices[support_point_index]] = mesh_displacement_values[q][d] - 
                                                                        old_mesh_displacement_values[q][d];
          }
      }

    distributed_vector.compress(VectorOperation::insert);
    mesh_displacement_increments = distributed_vector;
  }


  template <int dim>
  void MeshDeformationHandler<dim>::declare_parameters (ParameterHandler &prm)
  {
    prm.enter_subsection ("Mesh deformation");
    {
      prm.declare_entry ("Additional tangential mesh velocity boundary indicators", "",
                         Patterns::List(Patterns::Anything()),
                         "A comma separated list of names denoting those boundaries "
                         "where there the mesh is allowed to move tangential to the "
                         "boundary. All tangential mesh movements along "
                         "those boundaries that have tangential material velocity "
                         "boundary conditions are allowed by default, this parameters "
                         "allows to generate mesh movements along other boundaries that are "
                         "open, or have prescribed material velocities or tractions."
                         "\n\n"
                         "The names of the boundaries listed here can either be "
                         "numbers (in which case they correspond to the numerical "
                         "boundary indicators assigned by the geometry object), or they "
                         "can correspond to any of the symbolic names the geometry object "
                         "may have provided for each part of the boundary. You may want "
                         "to compare this with the documentation of the geometry model you "
                         "use in your model.");

      prm.declare_entry ("Free surface boundary indicators", "",
                         Patterns::List(Patterns::Anything()),
                         "A comma separated list of names denoting those boundaries "
                         "considered as free surfaces."
                         "\n\n"
                         "The names of the boundaries listed here can either be "
                         "numbers (in which case they correspond to the numerical "
                         "boundary indicators assigned by the geometry object), or they "
                         "can correspond to any of the symbolic names the geometry object "
                         "may have provided for each part of the boundary. You may want "
                         "to compare this with the documentation of the geometry model you "
                         "use in your model. ");

      prm.enter_subsection("Free surface");
      {
        prm.declare_entry("Free surface stabilization theta", "0.5",
                          Patterns::Double(0, 1),
                          "");
      }
      prm.leave_subsection();
    }
    prm.leave_subsection ();

    // Let the free surface model and mesh smoothing model
    // declare their own parameters.
    MeshDeformation::FreeSurface::declare_parameters<dim>(prm);
    MeshDeformation::MeshSmoothing::declare_parameters<dim>(prm);
  }


  template <int dim>
  void MeshDeformationHandler<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection ("Mesh deformation");
    {
      // Create the map of free surface boundary indicators
      try
      {
        const std::vector<types::boundary_id> x_free_surface_boundary_indicators
          = this->get_geometry_model().translate_symbolic_boundary_names_to_ids(
              Utilities::split_string_list(prm.get("Free surface boundary indicators")));

        free_surface_boundary_indicators.insert(x_free_surface_boundary_indicators.begin(),
                                                x_free_surface_boundary_indicators.end());
      }
      catch (const std::string &error)
      {
        AssertThrow(false, ExcMessage ("While parsing the entry <Mesh deformation/Free surface "
                                       "boundary indicators>, there was an error. Specifically, "
                                       "the conversion function complained as follows:\n\n"
                                       + error));
      }

      // Create the list of tangential mesh movement boundary indicators
      try
      {
        const std::vector<types::boundary_id> x_additional_tangential_mesh_boundary_indicators
          = this->get_geometry_model().translate_symbolic_boundary_names_to_ids(
              Utilities::split_string_list(prm.get("Additional tangential mesh velocity boundary indicators")));

        tangential_mesh_deformation_boundary_indicators.insert(x_additional_tangential_mesh_boundary_indicators.begin(),
                                                               x_additional_tangential_mesh_boundary_indicators.end());
      }
      catch (const std::string &error)
      {
        AssertThrow (false, ExcMessage ("While parsing the entry <Mesh deformation/Additional tangential "
                                        "mesh velocity boundary indicators>, there was an error. Specifically, "
                                        "the conversion function complained as follows:\n\n"
                                        + error));
      }

      // Tangential mesh boundaries and free surface boundaries should have no
      // intersection.
      std::vector<types::boundary_id> intersection(
        this->get_geometry_model().get_used_boundary_indicators().size());
      auto p = std::set_intersection(free_surface_boundary_indicators.begin(),
                                     free_surface_boundary_indicators.end(),
                                     tangential_mesh_deformation_boundary_indicators.begin(),
                                     tangential_mesh_deformation_boundary_indicators.end(),
                                     intersection.begin());
      AssertThrow(p == intersection.begin(),
                  ExcMessage("There is at least one boundary indicator that exists in both the "
                             "free surface boundary list and the additional tangential mesh "
                             "velocity boundary list!"));

      // Boundaries with tangential velocity boundary conditions are implicitly
      // treated as tangential mesh boundaries. Also, these boundaries should not
      // be treated as free surfaces.
      for (const auto &boundary_id : this->get_boundary_velocity_manager().get_tangential_boundary_velocity_indicators())
      {
        tangential_mesh_deformation_boundary_indicators.insert(boundary_id);
        free_surface_boundary_indicators.erase(boundary_id);
      }

      // Create the list of zero mesh movement (fixed) boundary indicators, these are
      // all boundaries, which are not free surface or tangential.
      zero_mesh_deformation_boundary_indicators = this->get_geometry_model().get_used_boundary_indicators();

      for (const auto &boundary_id : free_surface_boundary_indicators)
        zero_mesh_deformation_boundary_indicators.erase(boundary_id);

      for (const auto &boundary_id : tangential_mesh_deformation_boundary_indicators)
        zero_mesh_deformation_boundary_indicators.erase(boundary_id);

      prm.enter_subsection("Free surface");
      {
        free_surface_theta = prm.get_double("Free surface stabilization theta");
      }
      prm.leave_subsection();
    }
    prm.leave_subsection ();

    // Create the free surface model and the mesh smoothing model, and
    // let them parse their own parameters.
    free_surface_model.reset(MeshDeformation::FreeSurface::create_free_surface_model<dim>(prm));
    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(free_surface_model.get()))
      sim->initialize_simulator(this->get_simulator());
    free_surface_model->parse_parameters(prm);
    free_surface_model->initialize();

    mesh_smoothing_model.reset(MeshDeformation::MeshSmoothing::create_mesh_smoothing_model<dim>(prm));
    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(mesh_smoothing_model.get()))
      sim->initialize_simulator(this->get_simulator());
    mesh_smoothing_model->parse_parameters(prm);
    mesh_smoothing_model->initialize();

    // Initialize MeshDeformationHandler itself.
    mesh_deformation_dof_handler.reinit(this->get_triangulation());

    if (!Plugins::plugin_type_matches<InitialTopography::ZeroTopography<dim>>(this->get_initial_topography()))
      include_initial_topography = true;
  }


  template <int dim>
  const std::set<types::boundary_id> &
  MeshDeformationHandler<dim>::get_free_surface_boundary_indicators() const
  {
    return free_surface_boundary_indicators;
  }


  template <int dim>
  const TrilinosWrappers::MPI::BlockVector &
  MeshDeformationHandler<dim>::get_mesh_displacement_increments() const
  {
    return mesh_displacement_increments;
  }
}


// explicit instantiations
namespace elaspect
{
#define INSTANTIATE(dim) \
  template class MeshDeformationHandler<dim>; 

  ELASPECT_INSTANTIATE(INSTANTIATE)
}
