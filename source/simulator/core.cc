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
#include <elaspect/simulator/assemblers/interface.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/grid_refinement.h>

namespace elaspect
{
  template <int dim>
  Simulator<dim>::Simulator(const MPI_Comm mpi_communicator_,
                            ParameterHandler &prm)
    :
    assemblers(std::make_unique<Assemblers::Manager<dim> >()),
    parameters(prm, mpi_communicator_),
    introspection(construct_default_variables(parameters), parameters),
    mpi_communicator(Utilities::MPI::duplicate_communicator(mpi_communicator_)),
    iostream_tee_device(std::cout, log_file_stream),
    iostream_tee_stream(iostream_tee_device),
    pcout(iostream_tee_stream,
          (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
    statistics_last_write_size(0),
    statistics_last_hash(0),
    computing_timer(mpi_communicator,
                    pcout,
                    TimerOutput::never,
                    TimerOutput::wall_times),

    geometry_model(GeometryModel::create_geometry_model<dim>(prm)),
    initial_topography(InitialTopography::create_initial_topography_model<dim>(prm)),
    gravity_model(GravityModel::create_gravity_model<dim>(prm)),
    boundary_heat_flux(BoundaryHeatFlux::create_boundary_heat_flux<dim>(prm)),

    timestep_number(0),
    pre_refinement_step(0),
    nonlinear_iteration(0),

    triangulation(mpi_communicator,
                  typename Triangulation<dim>::MeshSmoothing
                  (
                    Triangulation<dim>::limit_level_difference_at_vertices |
                    (Triangulation<dim>::eliminate_unrefined_islands |
                     Triangulation<dim>::eliminate_refined_inner_islands |
                     Triangulation<dim>::do_not_produce_unrefined_islands)
                  )
                  ,
                  parallel::distributed::Triangulation<dim>::mesh_reconstruction_after_repartitioning),

    finite_element(introspection.get_fes(), introspection.get_multiplicities()),
    dof_handler(triangulation),
    quadrature_formula(parameters.n_gaussian_points),

    qpd_handler(triangulation, quadrature_formula.size(), introspection.n_qpd_components),

    rebuild_sparsity_and_matrices(true),
    reassemble_tangent_matrix(true)
  {
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      // only open the log file on processor 0, the other processors won't be
      // writing into the stream anyway
      log_file_stream.open((parameters.output_directory + "log.txt").c_str(), std::ios_base::out);

      // we already printed the header to the screen, so here we just dump it
      // into the log file.
      print_elaspect_header(log_file_stream);
    }

    // now that we have output set up, we can start timer sections
    TimerOutput::Scope timer (computing_timer, "Initialization");

    // if any plugin wants access to the Simulator by deriving from SimulatorAccess, initialize it and
    // call the initialize() functions immediately after.
    //
    // up front, we can not know whether a plugin derives from
    // SimulatorAccess. all we have is a pointer to the base class of
    // each plugin type (the 'Interface' class in the namespace
    // corresponding to each plugin type), but this base class is not
    // derived from SimulatorAccess. in order to find out whether a
    // concrete plugin derives from this base (interface) class AND
    // the SimulatorAccess class via multiple inheritance, we need to
    // do a sideways dynamic_cast to this putative sibling of the
    // interface class, and investigate if the dynamic_cast
    // succeeds. if it succeeds, the dynamic_cast returns a non-nullptr
    // result, and we can test this in an if-statement. there is a nice
    // idiom whereby we can write
    //    if (SiblingClass *ptr = dynamic_cast<SiblingClass*>(ptr_to_base))
    //      ptr->do_something()
    // where we declare a variable *inside* the 'if' condition, and only
    // enter the code block guarded by the 'if' in case the so-declared
    // variable evaluates to something non-zero, which here means that
    // the dynamic_cast succeeded and returned the address of the sibling
    // object.
    //
    // we also need to let all models parse their parameters. this is done *after* setting
    // up their SimulatorAccess base class so that they can query, for example, the
    // geometry model's description of symbolic names for boundary parts. note that
    // the geometry model is the only model whose run time parameters are already read
    // at the time it is created
    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(geometry_model.get()))
      sim->initialize_simulator(*this);
    geometry_model->initialize();

    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(initial_topography.get()))
      sim->initialize_simulator(*this);
    initial_topography->initialize();

    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(gravity_model.get()))
      sim->initialize_simulator(*this);
    gravity_model->parse_parameters(prm);
    gravity_model->initialize();

    if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(boundary_heat_flux.get()))
      sim->initialize_simulator(*this);

    boundary_heat_flux->parse_parameters(prm);
    boundary_heat_flux->initialize();

    material_handler.initialize_simulator(*this);
    material_handler.parse_parameters(prm);

    heating_model_manager.initialize_simulator(*this);
    heating_model_manager.parse_parameters(prm);

    initial_composition_manager.initialize_simulator(*this);
    initial_composition_manager.parse_parameters(prm);

    initial_temperature_manager.initialize_simulator(*this);
    initial_temperature_manager.parse_parameters(prm);

    boundary_composition_manager.initialize_simulator(*this);
    boundary_composition_manager.parse_parameters(prm);

    boundary_temperature_manager.initialize_simulator(*this);
    boundary_temperature_manager.parse_parameters(prm);

    boundary_velocity_manager.initialize_simulator(*this);
    boundary_velocity_manager.parse_parameters(prm);

    postprocess_manager.initialize_simulator(*this);
    postprocess_manager.parse_parameters(prm);

    mesh_refinement_manager.initialize_simulator(*this);
    mesh_refinement_manager.parse_parameters(prm);

    time_stepping_manager.initialize_simulator(*this);
    time_stepping_manager.parse_parameters(prm);

    geometry_model->create_coarse_mesh(triangulation);

    parameters.parse_geometry_dependent_parameters (prm, *geometry_model);

    mesh_deformation_handler.initialize_simulator(*this);
    mesh_deformation_handler.parse_parameters(prm);

    for (const auto &p : parameters.prescribed_traction_boundary_indicators)
    {
      BoundaryTraction::Interface<dim> *bv
        = BoundaryTraction::create_boundary_traction<dim>
          (p.second.second);
      boundary_traction[p.first].reset (bv);
      if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(bv))
        sim->initialize_simulator(*this);
      bv->parse_parameters (prm);
      bv->initialize ();
    }

    // make sure that we don't have to fill every column of the statistics
    // object in each time step.
    statistics.set_auto_fill_mode(true);

    set_assemblers();
  }


  template <int dim>
  Simulator<dim>::~Simulator ()
  {
    // wait if there is a thread that's still writing the statistics
    // object (set from the output_statistics() function)
    if (output_statistics_thread.joinable())
      output_statistics_thread.join();

    // If an exception is being thrown (for example due to AssertThrow()), we
    // might end up here with currently active timing sections. The destructor
    // of TimerOutput does MPI communication, which can lead to deadlocks,
    // hangs, or confusing MPI error messages. To avoid this, we can call
    // reset() to remove all open sections. In a normal run, we won't have any
    // active sessions, so this won't hurt to do:
    computing_timer.reset();
  }


  template <int dim>
  void Simulator<dim>::setup_introspection ()
  {
    introspection.system_dofs_per_block = DoFTools::count_dofs_per_fe_block (dof_handler,
                                                                             introspection.get_components_to_blocks());

    IndexSet system_index_set = dof_handler.locally_owned_dofs();
    introspection.index_sets.system_partitioning = system_index_set.split_by_block(introspection.system_dofs_per_block);
    
    DoFTools::extract_locally_relevant_dofs (dof_handler,
                                             introspection.index_sets.system_relevant_set);
    introspection.index_sets.system_relevant_partitioning =
      introspection.index_sets.system_relevant_set.split_by_block(introspection.system_dofs_per_block);
  }


  template <int dim>
  void Simulator<dim>::setup_system_matrix ()
  {
    system_matrix.clear();

    const typename Introspection<dim>::ComponentIndices &x
      = introspection.component_indices;

    // Set up the coupling table.
    Table<2,DoFTools::Coupling> coupling (introspection.n_components,
                                          introspection.n_components);
    coupling.fill (DoFTools::none);

    for (unsigned int c = 0; c < dim; ++c)
      for (unsigned int d = 0; d < dim; ++d)
        coupling[x.displacement[c]][x.displacement[d]] = DoFTools::always;

    if (parameters.include_heat_transport)
      coupling[x.temperature][x.temperature] = DoFTools::always;

    if (x.compositional_fields.size() > 0)
      coupling[x.compositional_fields[0]][x.compositional_fields[0]] = DoFTools::always;

    coupling[x.stress[0]][x.stress[0]] = DoFTools::always;

    if (parameters.constitutive_relation & ConstitutiveRelation::pore_fluid)
    {
      const unsigned int fp_component_index = 
        introspection.variable("fluid pressure").first_component_index;
      for (unsigned int d = 0; d < dim; ++d)
      {
        coupling[x.displacement[d]][fp_component_index] = DoFTools::always;
        coupling[fp_component_index][x.displacement[d]] = DoFTools::always;
      }

      coupling[fp_component_index][fp_component_index] = DoFTools::always;
    }

    TrilinosWrappers::BlockSparsityPattern
    sp (introspection.index_sets.system_partitioning,
        introspection.index_sets.system_partitioning,
        introspection.index_sets.system_relevant_partitioning,
        mpi_communicator);

    // Set up the face coupling table if ALE method is applied.
    if (parameters.use_ALE_method)
    {
      Table<2,DoFTools::Coupling> face_coupling(introspection.n_components,
                                                introspection.n_components);
      face_coupling.fill(DoFTools::none);

      if (parameters.include_heat_transport)
        face_coupling[x.temperature][x.temperature] = DoFTools::always;

      if (x.compositional_fields.size() > 0)
        face_coupling[x.compositional_fields[0]][x.compositional_fields[0]] = DoFTools::always;

      face_coupling[x.stress[0]][x.stress[0]] = DoFTools::always;

      DoFTools::make_flux_sparsity_pattern(dof_handler, sp,
                                           current_constraints, false,
                                           coupling,
                                           face_coupling,
                                           Utilities::MPI::
                                           this_mpi_process(mpi_communicator));
    }
    else
      DoFTools::make_sparsity_pattern (dof_handler,
                                       coupling, sp, 
                                       current_constraints, false,
                                       Utilities::MPI::
                                       this_mpi_process(mpi_communicator));

    sp.compress();

    // We may only allocate some of the matrix blocks, but the sparsity pattern
    // will still create entries for hanging nodes and boundary conditions.
    // These are unnecessary and are removed here.
    for (unsigned int i = 0; i < introspection.n_components; ++i)
      if (coupling[i][i] == DoFTools::none)
      {
        const unsigned int block = introspection.get_components_to_blocks()[i];

        sp.block(block,block).reinit(sp.block(block,block).locally_owned_range_indices(),
                                     sp.block(block,block).locally_owned_domain_indices());
        sp.block(block,block).compress();
      }

    system_matrix.reinit(sp);
  }


  template <int dim>
  void Simulator<dim>::setup_dofs ()
  {
    TimerOutput::Scope timer (computing_timer, "Setup dof systems");

    dof_handler.distribute_dofs(finite_element);

    // Renumber the DoFs hierarchical so that we get the
    // same numbering if we resume the computation. This
    // is because the numbering depends on the order the
    // cells are created.
    DoFRenumbering::hierarchical (dof_handler);
    DoFRenumbering::component_wise (dof_handler,
                                    introspection.get_components_to_blocks());

    // set up the introspection object that stores all sorts of
    // information about components of the finite element, component
    // masks, etc
    setup_introspection();

    // print dof numbers. Do so with 1000s separator since they are frequently
    // large
    {
      std::locale s = pcout.get_stream().getloc();
      // Creating std::locale with an empty string previously caused problems
      // on some platforms, so the functionality to catch the exception and ignore
      // is kept here, even though explicitly setting a facet should always work.
      try
      {
        // Imbue the stream with a locale that does the right thing. The
        // locale is responsible for later deleting the object pointed
        // to by the last argument (the "facet"), see
        // https://en.cppreference.com/w/cpp/locale/locale/locale
        pcout.get_stream().imbue(std::locale(std::locale(),
                                             new elaspect::Utilities::ThousandSep));
      }
      catch (const std::runtime_error &e)
      {
        // If the locale doesn't work, just give up
      }

      pcout << "Number of active cells: "
            << triangulation.n_global_active_cells()
            << " (on "
            << triangulation.n_global_levels()
            << " levels)"
            << std::endl
            << "Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << " ("
            << introspection.system_dofs_per_block[0];
      for (unsigned int b=1; b<introspection.system_dofs_per_block.size(); ++b)
        pcout << '+' << introspection.system_dofs_per_block[b];
      pcout <<')'
            << std::endl
            << std::endl;

      pcout.get_stream().imbue(s);
    }

    mapping = std::make_unique<MappingQ1Eulerian<dim, TrilinosWrappers::MPI::Vector>> (
      mesh_deformation_handler.mesh_deformation_dof_handler, 
      mesh_deformation_handler.mesh_displacements);

    mesh_deformation_handler.setup_dofs();
    pcout << std::endl;

    // Reconstruct the constraint-matrix:
    constraints.clear();
    constraints.reinit(introspection.index_sets.system_relevant_set);

    system_rhs.reinit (introspection.index_sets.system_partitioning,
                       mpi_communicator);
    solution.reinit (introspection.index_sets.system_partitioning,
                     introspection.index_sets.system_relevant_partitioning,
                     mpi_communicator);
    old_solution.reinit (introspection.index_sets.system_partitioning,
                         introspection.index_sets.system_relevant_partitioning,
                         mpi_communicator);

    //  Make hanging node constraints:
    DoFTools::make_hanging_node_constraints (dof_handler,
                                             constraints);

    compute_initial_displacement_boundary_constraints(constraints);
    constraints.close();

    qpd_handler.reinit();
  } 


  template <int dim>
  void 
  Simulator<dim>::
  compute_initial_displacement_boundary_constraints (AffineConstraints<double> &constraints)
  {
    // This needs to happen after the periodic constraints are added:
    setup_nullspace_constraints(constraints);

    // do the interpolation for zero displacement
    for (const auto p : boundary_velocity_manager.get_zero_boundary_velocity_indicators())
      VectorTools::interpolate_boundary_values (*mapping,
                                                dof_handler,
                                                p,
                                                Functions::ZeroFunction<dim>(introspection.n_components),
                                                constraints,
                                                introspection.component_masks.displacement);

    // do the same for no-normal-flux boundaries
    VectorTools::compute_no_normal_flux_constraints (dof_handler,
                                                     /* first_vector_component= */
                                                     introspection.component_indices.displacement[0],
                                                     boundary_velocity_manager.get_tangential_boundary_velocity_indicators(),
                                                     constraints,
                                                     *mapping);
  }


  template <int dim>
  void Simulator<dim>::compute_current_constraints ()
  {
    // We put the constraints we compute into a separate AffineConstraints<double> so we can check
    // if the set of constraints has changed. If it did, we need to update the sparsity patterns.
    AffineConstraints<double> new_current_constraints;
    new_current_constraints.clear();
    new_current_constraints.reinit (introspection.index_sets.system_relevant_set);
    new_current_constraints.merge (constraints);
    compute_current_displacement_boundary_constraints(new_current_constraints);

    // If there is a fixed boundary temperature or heat flux,
    // update the temperature boundary condition.
    boundary_temperature_manager.update();
    boundary_heat_flux->update();

    // Do the same for compositional fields.
    boundary_composition_manager.update();

    // If ALE method is applied, then both temperature and compositional
    // fields are discretized with discontinuous FE and boundary conditions
    // are not handled with AffineConstraints.
    if (!parameters.use_ALE_method)
    {
      // obtain the boundary indicators that belong to Dirichlet-type
      // temperature boundary conditions and interpolate the temperature
      // there
      for (const auto p : boundary_temperature_manager.get_fixed_temperature_boundary_indicators())
      {
        auto lambda = [&] (const Point<dim> &x) -> double
        {
          return boundary_temperature_manager.boundary_temperature(p, x);
        };

        VectorFunctionFromScalarFunctionObject<dim> vector_function_object(
          lambda,
          introspection.component_masks.temperature.first_selected_component(),
          introspection.n_components);

        VectorTools::interpolate_boundary_values (*mapping,
                                                  dof_handler,
                                                  p,
                                                  vector_function_object,
                                                  new_current_constraints,
                                                  introspection.component_masks.temperature);
      }
    }

    new_current_constraints.close();

    // Now check if the current_constraints we just computed changed from before. 
    // If the mesh got refined and the size of the linear system changed, the old
    // and new constraint matrices will have different entries, and we can not 
    // easily compare them.
    const bool mesh_has_changed = (current_constraints.get_local_lines().size()
                                   != new_current_constraints.get_local_lines().size())
                                  ||
                                  (current_constraints.get_local_lines()
                                   != new_current_constraints.get_local_lines());
    // Figure out if any entry that was constrained before is now no longer constrained:
    bool constraints_changed = false;

    if (mesh_has_changed)
    {
      constraints_changed = true;
    }
    else
    {
      // the mesh has not changed on our machine, so compare the constraints:
      for (auto &row: current_constraints.get_lines())
      {
        if (!new_current_constraints.is_constrained(row.index))
        {
          constraints_changed = true;
          break;
        }
      }
    }

    // If at least one processor has different constraints, force rebuilding the matrices:
    const bool any_constrained_dofs_set_changed = Utilities::MPI::sum(constraints_changed ? 1 : 0,
                                                                      mpi_communicator) > 0;
    if (any_constrained_dofs_set_changed)
      rebuild_sparsity_and_matrices = true;

    current_constraints.copy_from(new_current_constraints);
  }


  template <int dim>
  void 
  Simulator<dim>::
  compute_current_displacement_boundary_constraints (AffineConstraints<double> &constraints)
  {
    boundary_velocity_manager.update();
    for (const auto &p : boundary_velocity_manager.get_active_boundary_velocity_names())
    {
      Utilities::VectorFunctionFromVelocityFunctionObject<dim> disp
      (introspection.n_components, 
       [&] (const Point<dim> &x) -> Tensor<1,dim>
      {
        return boundary_velocity_manager.boundary_velocity(p.first, x) * time_step;
      });

      std::vector<bool> mask(introspection.component_masks.displacement.size(), false);
      const std::string &comp = p.second.first;

      if (comp.length() > 0)
      {
        for (std::string::const_iterator direction = comp.begin(); direction != comp.end(); ++direction)
        {
          switch (*direction)
          {
            case 'x':
              mask[introspection.component_indices.displacement[0]] = true;
              break;
            case 'y':
              mask[introspection.component_indices.displacement[1]] = true;
              break;
            case 'z':
              // we must be in 3d, or 'z' should never have gotten through
              Assert (dim==3, ExcInternalError());
              if (dim==3)
                mask[introspection.component_indices.displacement[dim-1]] = true;
              break;
            default:
              Assert (false, ExcInternalError());
          }
        }
      }
      else
      {
        // no mask given -- take all displacements
        for (unsigned int i=0; i<introspection.component_masks.displacement.size(); ++i)
          mask[i]=introspection.component_masks.displacement[i];
      }

      if (!(parameters.constitutive_relation & ConstitutiveRelation::plasticity)
          || (nonlinear_iteration == 0))
      {
        VectorTools::interpolate_boundary_values (*mapping,
                                                  dof_handler,
                                                  p.first,
                                                  disp,
                                                  constraints,
                                                  mask);
      }
      else
        VectorTools::interpolate_boundary_values (*mapping,
                                                  dof_handler,
                                                  p.first,
                                                  Functions::ZeroFunction<dim>(introspection.n_components),
                                                  constraints,
                                                  mask);
    }
  }


  template <int dim>
  void
  Simulator<dim>::start_timestep ()
  {
    nonlinear_iteration = 0;

    // produce some output for the screen to show where we are
    const char *unit = (parameters.convert_to_years ? "years" : "seconds");
    const double multiplier = (parameters.convert_to_years ? 1./year_in_seconds : 1.0);

    pcout << "*** Timestep " << timestep_number
          << ":  t=" << (time * multiplier) << ' ' << unit
          << ", dt=" << (time_step * multiplier) << ' ' << unit
          << std::endl;

    // interpolate the current boundary displacements. copy constraints
    // into current_constraints and then add to current_constraints
    compute_current_constraints();

    // If needed, construct sparsity patterns and matrices with the current
    // constraints. Of course we need to force assembly too.
    if (rebuild_sparsity_and_matrices)
    {
      TimerOutput::Scope timer (computing_timer, "Setup matrices");

      setup_system_matrix();
      rebuild_sparsity_and_matrices = false;
      reassemble_tangent_matrix = true;
    }

    gravity_model->update();
    material_handler.update();
    mesh_refinement_manager.update();
    mesh_deformation_handler.update();
  }


  template <int dim>
  void Simulator<dim>::postprocess ()
  {
    TimerOutput::Scope timer (computing_timer, "postprocessing");
    pcout << "   Postprocessing:" << std::endl;

    // run all the postprocessing routines and then write
    // the current state of the statistics table to a file
    std::list<std::pair<std::string,std::string> >
    output_list = postprocess_manager.execute (statistics);

    // if we are on processor zero, print to screen
    // whatever the postprocessors have generated
    if (Utilities::MPI::this_mpi_process(mpi_communicator)==0)
    {
      // determine the width of the first column of text so that
      // everything gets nicely aligned; then output everything
      {
        unsigned int width = 0;
        for (const auto &p : output_list)
          width = std::max<unsigned int> (width, p.first.size());

        for (const auto &p : output_list)
          pcout << "     "
                << std::left
                << std::setw(width)
                << p.first
                << " "
                << p.second
                << std::endl;
      }

      pcout << std::endl;
    }

    // finally, write the entire set of current results to disk
    output_statistics();
  }


  template <int dim>
  void Simulator<dim>::refine_mesh (const unsigned int max_grid_level)
  {
    parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::BlockVector>
    system_trans(dof_handler);

    parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector>
    md_system_trans(mesh_deformation_handler.mesh_deformation_dof_handler);

    {
      TimerOutput::Scope timer (computing_timer, "Refine mesh structure, part 1");

      Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
      mesh_refinement_manager.execute (estimated_error_per_cell);

      if (parameters.adapt_by_fraction_of_cells)
        parallel::distributed::GridRefinement::
        refine_and_coarsen_fixed_number (triangulation,
                                         estimated_error_per_cell,
                                         parameters.refinement_fraction,
                                         parameters.coarsening_fraction);
      else
        parallel::distributed::GridRefinement::
        refine_and_coarsen_fixed_fraction (triangulation,
                                           estimated_error_per_cell,
                                           parameters.refinement_fraction,
                                           parameters.coarsening_fraction);

      mesh_refinement_manager.tag_additional_cells();

      // clear refinement flags if parameter.refinement_fraction=0.0
      if (parameters.refinement_fraction == 0.0)
      {
        for (const auto &cell : triangulation.active_cell_iterators())
          cell->clear_refine_flag();
      }
      else
      {
        // limit maximum refinement level
        if (triangulation.n_levels() > max_grid_level)
          for (typename Triangulation<dim>::active_cell_iterator
               cell = triangulation.begin_active(max_grid_level);
               cell != triangulation.end(); ++cell)
            cell->clear_refine_flag();
      }

      // clear coarsening flags if parameter.coarsening_fraction=0.0
      if (parameters.coarsening_fraction == 0.0)
      {
        for (const auto &cell : triangulation.active_cell_iterators())
          cell->clear_coarsen_flag();
      }
      else
      {
        // limit miminum refinement level
        for (typename Triangulation<dim>::active_cell_iterator
             cell = triangulation.begin_active(0);
             cell != triangulation.end_active(parameters.min_grid_level); ++cell)
          cell->clear_coarsen_flag();
      }

      std::vector<const TrilinosWrappers::MPI::BlockVector *> x_system(1);
      x_system[0] = &solution;
      if (parameters.use_ALE_method)
        x_system.push_back(&mesh_deformation_handler.mesh_displacement_increments);

      std::vector<const TrilinosWrappers::MPI::Vector *> x_md_system(2);
      x_md_system[0] = &mesh_deformation_handler.mesh_displacements;
      x_md_system[1] = &mesh_deformation_handler.initial_topography;

      {
        // Communicate refinement flags on ghost cells from the owner of the
        // cell. This is necessary to get consistent refinement, as mesh
        // smoothing would undo some of the requested coarsening/refinement.
        auto pack
        = [] (const typename DoFHandler<dim>::active_cell_iterator &cell) -> std::uint8_t
        {
          if (cell->refine_flag_set())
            return 1;
          if (cell->coarsen_flag_set())
            return 2;
          return 0;
        };

        auto unpack
        = [] (const typename DoFHandler<dim>::active_cell_iterator &cell, const std::uint8_t &flag) -> void
        {
          cell->clear_coarsen_flag();
          cell->clear_refine_flag();
          if (flag==1)
            cell->set_refine_flag();
          else if (flag==2)
            cell->set_coarsen_flag();
        };

        GridTools::exchange_cell_data_to_ghosts<std::uint8_t, DoFHandler<dim>> (dof_handler, pack, unpack);
      }

      triangulation.prepare_coarsening_and_refinement();
      system_trans.prepare_for_coarsening_and_refinement(x_system);
      md_system_trans.prepare_for_coarsening_and_refinement(x_md_system);

      triangulation.execute_coarsening_and_refinement();
    } // leave the timed section

    setup_dofs ();

    {
      TimerOutput::Scope timer (computing_timer, "Refine mesh structure, part 2");

      TrilinosWrappers::MPI::BlockVector distributed_system;
      TrilinosWrappers::MPI::BlockVector distributed_mesh_displacement_increments;
      distributed_system.reinit(introspection.index_sets.system_partitioning, mpi_communicator);
      if (parameters.use_ALE_method)
        distributed_mesh_displacement_increments.reinit(
            introspection.index_sets.system_partitioning, mpi_communicator);

      std::vector<TrilinosWrappers::MPI::BlockVector *> system_tmp(1);
      system_tmp[0] = &distributed_system;
      if (parameters.use_ALE_method)
        system_tmp.push_back(&distributed_mesh_displacement_increments);

      // transfer the data previously stored into the vectors indexed by
      // system_tmp. then ensure that the interpolated solution satisfies
      // hanging node constraints
      //
      // note that the 'constraints' variable contains hanging node constraints
      // and constraints from periodic boundary conditions (as well as from
      // zero and tangential velocity boundary conditions), but not from
      // non-homogeneous boundary conditions. the latter are added not in setup_dofs(),
      // which we call above, but added to 'current_constraints' in start_timestep(),
      // which we do not want to call here.
      //
      // however, what we have should be sufficient: we have everything that
      // is necessary to make the solution vectors *conforming* on the current
      // mesh.
      system_trans.interpolate (system_tmp);

      constraints.distribute (distributed_system);
      solution = distributed_system;

      // refresh the quadrature point data with the solution vector
      refresh_quadrature_point_data();

      if (parameters.use_ALE_method)
      {
        constraints.distribute (distributed_mesh_displacement_increments);
        mesh_deformation_handler.mesh_displacement_increments = 
          distributed_mesh_displacement_increments;
      }

      // do the same as above, but for the mesh deformation vectors
      TrilinosWrappers::MPI::Vector distributed_mesh_displacements;
      TrilinosWrappers::MPI::Vector distributed_initial_topography;
      distributed_mesh_displacements.reinit (mesh_deformation_handler.mesh_locally_owned,
                                             mpi_communicator);
      distributed_initial_topography.reinit (mesh_deformation_handler.mesh_locally_owned,
                                             mpi_communicator);

      std::vector<TrilinosWrappers::MPI::Vector *> md_system_tmp(2);
      md_system_tmp[0] = &distributed_mesh_displacements;
      md_system_tmp[1] = &distributed_initial_topography;

      md_system_trans.interpolate (md_system_tmp);

      mesh_deformation_handler.mesh_vertex_constraints.distribute (distributed_mesh_displacements);
      mesh_deformation_handler.mesh_displacements = distributed_mesh_displacements;

      mesh_deformation_handler.mesh_vertex_constraints.distribute (distributed_initial_topography);
      mesh_deformation_handler.initial_topography = distributed_initial_topography;

      // calculate global volume after refining mesh
      global_volume = GridTools::volume(triangulation, *mapping);
    }
  }


  template <int dim>
  void Simulator<dim>::run ()
  {
    unsigned int max_refinement_level = parameters.initial_global_refinement +
                                        parameters.initial_adaptive_refinement;

    pre_refinement_step = 0;

    time = parameters.start_time;

    triangulation.refine_global(parameters.initial_global_refinement);

    setup_dofs();

    // calculate global volume after refining mesh
    global_volume = GridTools::volume(triangulation, *mapping);

start_time_iteration:

    TimerOutput::Scope timer (computing_timer, "Setup initial conditions");

    timestep_number = 0;

    // Add topography after all initial refinement is completed.
    if (pre_refinement_step == parameters.initial_adaptive_refinement)
      signals.pre_set_initial_state(triangulation);

    initialize_temperature_field();
    initialize_fluid_pressure_field();
    initialize_compositional_fields();

    time_stepping_manager.update();
    time_step = time_stepping_manager.get_next_time_step_size();

    do
    {
      if (! (parameters.skip_solvers_on_initial_refinement
             && pre_refinement_step < parameters.initial_adaptive_refinement))
      {
        start_timestep();

        do
        {
          // Solve the mechanical system and check the convergence.
          bool converged = assemble_and_solve_mechanical_system();

          if (!converged && parameters.enforce_convergence_for_mechanical_system)
          {
            // If the mechanical system failed to converge, check if we can
            // cut the time step size and try again.
            if (time_step > time_stepping_manager.get_minimum_time_step_size())
            {
              pcout << "*** The mechanical system failed to converge. "
                    << std::endl
                    << "*** We shall cut the time step by half and try again."
                    << std::endl
                    << std::endl;

              time_step *= 0.5;

              // We have to call start_timestep() again because some of the
              // plugins may depend on the time step size.
              start_timestep();
            }
            else
            {
              pcout << "ERROR: The mechanical system failed to converge with the "
                    << "smallest possible time step size. "
                    << std::endl;

              throw elaspect::QuietException();
            }
          }
          else
          {
            // If the mechanical system is converged or the convergence of 
            // the mechanical system is not essential, break the loop.
            break;
          }
        }
        while (true);
      }

      pcout << std::endl;

      // See if we have to start over with a new adaptive refinement cycle
      // at the beginning of the simulation.
      if (timestep_number == 0)
      {
        const bool initial_refinement_done = maybe_do_initial_refinement(max_refinement_level);
        if (initial_refinement_done)
          goto start_time_iteration;
      }

      // If we postprocess nonlinear iterations, this function is called within
      // assemble_and_solve_mechanical_system() in the individual solver schemes
      if (!parameters.run_postprocessors_on_nonlinear_iterations)
        postprocess();

      // Update the mesh deformation handler and execute mesh deformation.
      mesh_deformation_handler.execute();
      
      // calculate global volume after deforming mesh
      global_volume = GridTools::volume(triangulation, *mapping);

      // Update the data stored on quadrature points and project/advect them
      // onto/across grid nodes.
      update_quadrature_point_data();
      assemble_and_solve_qpd_system();
      refresh_quadrature_point_data();

      // Solve the thermo- and hydro-system.
      assemble_and_solve_thermo_system();
      assemble_and_solve_hydro_system();

      pcout << std::endl;

      // see if we need to refine the mesh
      maybe_refine_mesh(max_refinement_level);

      // see if we want to write a timing summary
      maybe_write_timing_output();

      // advance time
      time += time_step;
      ++timestep_number;

      if (time_stepping_manager.should_simulation_terminate_now())
        break;

      // update time step
      time_stepping_manager.update();
      time_step = time_stepping_manager.get_next_time_step_size();

      old_solution = solution;
    }
    while (true);

    // we disable automatic summary printing so that it won't happen when
    // throwing an exception. Therefore, we have to do this manually here:
    computing_timer.print_summary();
  }

}


// explicit instantiations
namespace elaspect
{
#define INSTANTIATE(dim) \
  template class Simulator<dim>;

  ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
