#include <elaspect/simulator.h>
#include <elaspect/simulator/assemblers/interface.h>
#include <elaspect/simulator/assemblers/thermo_system.h>
#include <elaspect/simulator/assemblers/hydro_system.h>
#include <elaspect/simulator/assemblers/mechanical_system.h>
#include <elaspect/simulator/assemblers/qpd_system.h>

#include <deal.II/base/work_stream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/physics/elasticity/standard_tensors.h>

namespace elaspect
{
  namespace
  {
    // This function initializes the simulator access for all assemblers
    // inside of the assemblers parameter. This is just a shortcut to save some
    // lines in set_assemblers(), where this operation appears multiple times.
    template <int dim, class AssemblerType>
    void
    initialize_simulator(const Simulator<dim> &simulator,
                         std::vector<std::unique_ptr<AssemblerType>> &assemblers)
    {
      for (unsigned int i=0; i<assemblers.size(); ++i)
        if (SimulatorAccess<dim> *p = dynamic_cast<SimulatorAccess<dim>*>(assemblers[i].get()))
          p->initialize_simulator(simulator);
    }
  }


  template <int dim>
  void Simulator<dim>::set_assemblers ()
  {
    // first let the manager delete all existing assemblers:
    assemblers->reset();

    assemblers->thermo_system.push_back(std::make_unique<Assemblers::ThermoConvectionDiffusion<dim>>());
    if (parameters.fixed_heat_flux_boundary_indicators.size() > 0)
      assemblers->thermo_system.push_back(std::make_unique<Assemblers::ThermoBoundaryFlux<dim>>());

    if (parameters.constitutive_relation & ConstitutiveRelation::pore_fluid)
      assemblers->hydro_system.push_back(std::make_unique<Assemblers::HydroUnsaturatedFlow<dim>>());

    assemblers->mechanical_system.push_back(std::make_unique<Assemblers::MechanicalQuasiStaticTerms<dim>>());

    // add the terms for traction boundary conditions
    if (!boundary_traction.empty())
      assemblers->mechanical_system_on_boundary_face.push_back(
        std::make_unique<Assemblers::MechanicalBoundaryTraction<dim>>());

    if (parameters.use_ALE_method)
    {
      assemblers->qpd_system.push_back(std::make_unique<Assemblers::QPDAdvection<dim>>());
      assemblers->qpd_system_on_boundary_face.push_back(
        std::make_unique<Assemblers::QPDAdvectionBoundaryFace<dim>>());
      assemblers->qpd_system_on_interior_face.push_back(
        std::make_unique<Assemblers::QPDAdvectionInteriorFace<dim>>());
    }
    else
      assemblers->qpd_system.push_back(std::make_unique<Assemblers::QPDProjection<dim>>());

    // allow other assemblers to add themselves or modify the existing ones by firing the signal
    signals.set_assemblers(*this, *assemblers);

    // ensure that all assembler objects have access to the SimulatorAccess
    initialize_simulator(*this, assemblers->mechanical_system);
    initialize_simulator(*this, assemblers->mechanical_system_on_boundary_face);
    initialize_simulator(*this, assemblers->qpd_system);
    initialize_simulator(*this, assemblers->qpd_system_on_boundary_face);
    initialize_simulator(*this, assemblers->qpd_system_on_interior_face);
  }


  template <int dim>
  void 
  Simulator<dim>::
  local_assemble_mechanical_system (const typename DoFHandler<dim>::active_cell_iterator &cell,
                               internal::Assembly::Scratch::MechanicalSystem<dim> &scratch,
                               internal::Assembly::CopyData::MechanicalSystem<dim> &data)
  {
    scratch.fe_values.reinit(cell);
    scratch.cell = cell;

    scratch.material_model_inputs.reinit (scratch.fe_values, 
                                          qpd_handler,
                                          introspection,
                                          solution);

    for (unsigned int i = 0; i < assemblers->mechanical_system.size(); ++i)
      assemblers->mechanical_system[i]->create_additional_material_model_outputs(scratch.material_model_outputs);

    material_handler.evaluate (scratch.material_model_inputs,
                               scratch.material_model_outputs);
    MaterialModel::MaterialAveraging::average (parameters.material_averaging,
                                               cell,
                                               scratch.fe_values.get_quadrature(),
                                               scratch.fe_values.get_mapping(),
                                               (parameters.constitutive_relation & ConstitutiveRelation::plasticity),
                                               scratch.material_model_outputs);

    // get all dof indices in the current cell, then extract those 
    // that correspond to the solution fields we are interested in.
    cell->get_dof_indices (scratch.local_dof_indices);

    const unsigned int u_dofs_per_cell = data.local_dof_indices.size();
    for (unsigned int i = 0, i_u = 0; i_u < u_dofs_per_cell; /*increment at end of loop*/)
    {
      if (finite_element.system_to_component_index(i).first < dim)
      {
        data.local_dof_indices[i_u] = scratch.local_dof_indices[i];
        ++i_u;
      }
      ++i;
    }

    data.local_matrix = 0;
    data.local_rhs = 0;

    // trigger the invocation of the various functions that actually do
    // all of the assembling
    for (unsigned int i = 0; i < assemblers->mechanical_system.size(); ++i)
      assemblers->mechanical_system[i]->execute(scratch, data);

    if (!assemblers->mechanical_system_on_boundary_face.empty())
    {
      // then also work on possible face terms.
      for (const unsigned int face_no : cell->face_indices())
        if (cell->at_boundary(face_no) && !cell->has_periodic_neighbor(face_no))
        {
          scratch.fe_face_values.reinit(cell, face_no);
          scratch.face_number = face_no;
          for (unsigned int i = 0; i < assemblers->mechanical_system_on_boundary_face.size(); ++i)
            assemblers->mechanical_system_on_boundary_face[i]->execute(scratch, data);
        }
    }
  }


  template <int dim>
  void 
  Simulator<dim>::
  copy_local_to_global_mechanical_system (const internal::Assembly::CopyData::MechanicalSystem<dim> &data)
  {
    current_constraints.distribute_local_to_global (data.local_matrix,
                                                    data.local_rhs,
                                                    data.local_dof_indices,
                                                    system_matrix,
                                                    system_rhs);
  }


  template <int dim>
  void Simulator<dim>::assemble_mechanical_system ()
  {
    TimerOutput::Scope timer (computing_timer,
                              "Assemble mechanical system");

    system_rhs    = 0;
    system_matrix = 0;

    auto worker = [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
                      internal::Assembly::Scratch::MechanicalSystem<dim> &scratch,
                      internal::Assembly::CopyData::MechanicalSystem<dim> &data)
    {
      this->local_assemble_mechanical_system (cell, scratch, data);
    };

    auto copier = [&](const internal::Assembly::CopyData::MechanicalSystem<dim> &data)
    {
      this->copy_local_to_global_mechanical_system (data);
    };

    const QGauss<dim-1> face_quadrature_formula (parameters.n_gaussian_points);

    const UpdateFlags update_flags = update_values |
                                     update_gradients |
                                     update_quadrature_points |
                                     update_JxW_values;

    // TODO: 
    const UpdateFlags face_update_flags = 
      (!assemblers->mechanical_system_on_boundary_face.empty()
       ? update_flags | update_normal_vectors
       : update_default);

    MaterialModel::FieldDependences::Dependence field_dependences = MaterialModel::FieldDependences::none;
    MaterialModel::MaterialProperties::Property requested_properties = MaterialModel::MaterialProperties::none;

    for (const auto &p : assemblers->mechanical_system)
    {
      field_dependences    |= p->get_field_dependences();
      requested_properties |= p->get_needed_material_properties();
    }
    for (const auto &p : assemblers->mechanical_system_on_boundary_face)
    {
      field_dependences   |= p->get_field_dependences();
      requested_properties |= p->get_needed_material_properties();
    }
    field_dependences |= material_handler.get_field_dependences_for_evaluation (requested_properties);

    const unsigned int u_dofs_per_cell = dim * finite_element.base_element(introspection.base_elements.displacement).dofs_per_cell;

    using CellFilter = FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>;

    WorkStream::
    run(CellFilter(IteratorFilters::LocallyOwnedCell(),
                   dof_handler.begin_active()),
        CellFilter(IteratorFilters::LocallyOwnedCell(),
                   dof_handler.end()),
        worker,
        copier,
        internal::Assembly::Scratch::MechanicalSystem<dim>(finite_element, 
                                                      *mapping,
                                                      quadrature_formula, 
                                                      face_quadrature_formula,
                                                      update_flags,
                                                      face_update_flags,
                                                      u_dofs_per_cell,
                                                      introspection.n_compositional_fields,
                                                      field_dependences,
                                                      requested_properties),
         internal::Assembly::CopyData::MechanicalSystem<dim> (u_dofs_per_cell));

    system_matrix.compress (VectorOperation::add);
    system_rhs.compress (VectorOperation::add);
  }


  template <int dim>
  void
  Simulator<dim>::
  local_assemble_thermo_system (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                 internal::Assembly::Scratch::ThermoSystem<dim> &scratch,
                                 internal::Assembly::CopyData::ThermoSystem<dim> &data)
  {
    scratch.fe_values.reinit(cell);
    scratch.cell = cell;

    // get all dof indices on the current cell, then extract those
    // that correspond to the temperature field
    cell->get_dof_indices (scratch.local_dof_indices);
    const unsigned int T_dofs_per_cell = data.local_dof_indices.size();
    for (unsigned int i = 0, i_T = 0; i_T < T_dofs_per_cell; /*increment at end of loop*/)
    {
      if (finite_element.system_to_component_index(i).first == introspection.component_indices.temperature)
      {
        data.local_dof_indices[i_T] = scratch.local_dof_indices[i];
        ++i_T;
      }
      ++i;
    }

    scratch.fe_values[introspection.extractors.temperature].get_function_values(old_solution,
                                                                                scratch.old_temperature_values);
    if (parameters.constitutive_relation & ConstitutiveRelation::pore_fluid)
    {
      const FEValuesExtractors::Scalar fp_extractor(introspection.variable("fluid pressure").first_component_index);
      scratch.fe_values[fp_extractor].get_function_gradients(old_solution, scratch.fluid_pressure_gradients);
    }

    if (parameters.use_ALE_method)
    {
      scratch.fe_values[introspection.extractors.displacement].get_function_values(solution,
                                                                                   scratch.displacement_increments);
      scratch.fe_values[introspection.extractors.displacement].get_function_values(mesh_deformation_handler.mesh_displacement_increments,
                                                                                   scratch.mesh_displacement_increments);
    }

    scratch.material_model_inputs.reinit (scratch.fe_values,
                                          qpd_handler,
                                          introspection,
                                          solution);

    for (unsigned int i = 0; i < assemblers->thermo_system.size(); ++i)
      assemblers->thermo_system[i]->create_additional_material_model_outputs(scratch.material_model_outputs);

    material_handler.evaluate (scratch.material_model_inputs,
                               scratch.material_model_outputs);
    MaterialModel::MaterialAveraging::average (parameters.material_averaging,
                                               cell,
                                               scratch.fe_values.get_quadrature(),
                                               scratch.fe_values.get_mapping(),
                                               false,
                                               scratch.material_model_outputs);

    heating_model_manager.evaluate(scratch.material_model_inputs,
                                   scratch.material_model_outputs,
                                   scratch.heating_model_outputs);

    data.local_matrix = 0;
    data.local_rhs = 0;
    
    for (unsigned int i = 0; i < assemblers->thermo_system.size(); ++i)
      assemblers->thermo_system[i]->execute(scratch, data);

    if (!assemblers->thermo_system_on_boundary_face.empty())
    {
      for (unsigned int face_number = 0; face_number < GeometryInfo<dim>::faces_per_cell; ++face_number)
      {
        const typename DoFHandler<dim>::face_iterator face = cell->face(face_number);
        if (face->at_boundary())
        {
          scratch.fe_face_values.reinit(cell, face_number);
          scratch.face_number = face_number;

          for (unsigned int i = 0; i < assemblers->thermo_system_on_boundary_face.size(); ++i)
            assemblers->thermo_system_on_boundary_face[i]->execute (scratch, data);
        }
      }
    }
  }


  template <int dim>
  void
  Simulator<dim>::
  copy_local_to_global_thermo_system (const internal::Assembly::CopyData::ThermoSystem<dim> &data)
  {
    current_constraints.distribute_local_to_global (data.local_matrix,
                                                    data.local_rhs,
                                                    data.local_dof_indices,
                                                    system_matrix,
                                                    system_rhs);
  }


  template <int dim>
  void Simulator<dim>::assemble_thermo_system ()
  {
    TimerOutput::Scope timer (computing_timer,
                              "Assemble thermo system");

    const unsigned int block_idx = introspection.block_indices.temperature;
    system_rhs.block (block_idx) = 0;
    system_matrix.block (block_idx, block_idx) = 0;

    auto worker = [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
                      internal::Assembly::Scratch::ThermoSystem<dim> &scratch,
                      internal::Assembly::CopyData::ThermoSystem<dim> &data)
    {
      this->local_assemble_thermo_system (cell, scratch, data);
    };

    auto copier = [&](const internal::Assembly::CopyData::ThermoSystem<dim> &data)
    {
      this->copy_local_to_global_thermo_system (data);
    };

    const QGauss<dim-1> face_quadrature_formula (parameters.n_gaussian_points);

    const UpdateFlags update_flags = update_values |
                                     update_gradients |
                                     update_quadrature_points |
                                     update_JxW_values;

    const UpdateFlags face_update_flags =
      (!assemblers->thermo_system_on_boundary_face.empty()
       ?
       update_values |
       update_quadrature_points |
       update_normal_vectors |
       update_JxW_values 
       :
       update_default);

    MaterialModel::FieldDependences::Dependence field_dependences = MaterialModel::FieldDependences::none;
    MaterialModel::MaterialProperties::Property requested_properties = MaterialModel::MaterialProperties::none;

    for (const auto &p : assemblers->thermo_system)
    {
      field_dependences    |= p->get_field_dependences();
      requested_properties |= p->get_needed_material_properties();
    }
    for (const auto &p : assemblers->thermo_system_on_boundary_face)
    {
      field_dependences    |= p->get_field_dependences();
      requested_properties |= p->get_needed_material_properties();
    }
    field_dependences |= material_handler.get_field_dependences_for_evaluation (requested_properties);
      
    const unsigned int T_dofs_per_cell = finite_element.base_element(introspection.base_elements.temperature).dofs_per_cell;

    using CellFilter = FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>;

    WorkStream::
    run(CellFilter(IteratorFilters::LocallyOwnedCell(),
                   dof_handler.begin_active()),
        CellFilter(IteratorFilters::LocallyOwnedCell(),
                   dof_handler.end()),
        worker,
        copier,
        internal::Assembly::Scratch::ThermoSystem<dim>(finite_element,
                                                        *mapping,
                                                        quadrature_formula,
                                                        face_quadrature_formula,
                                                        update_flags,
                                                        face_update_flags,
                                                        T_dofs_per_cell,
                                                        parameters.n_compositional_fields,
                                                        field_dependences,
                                                        requested_properties),
        internal::Assembly::CopyData::ThermoSystem<dim> (T_dofs_per_cell));

    system_matrix.compress (VectorOperation::add);
    system_rhs.compress (VectorOperation::add);
  }


  template <int dim>
  void
  Simulator<dim>::
  local_assemble_qpd_system(const typename DoFHandler<dim>::active_cell_iterator &cell,
                            internal::Assembly::Scratch::QPDSystem<dim> &scratch,
                            internal::Assembly::CopyData::QPDSystem<dim> &data)
  {
    scratch.fe_values.reinit(cell);
    scratch.cell = cell;

    const unsigned int n_fields = scratch.qpd_fields.size();
    const unsigned int field_dofs_per_cell = data.local_dof_indices[0].size();

    // get all dof indices in the current cell, then extract those
    // that correspond to the solution fields we are interested in.
    cell->get_dof_indices(scratch.local_dof_indices);
    for (unsigned int field = 0; field < scratch.qpd_fields.size(); ++field)
    {
      const unsigned int field_component = scratch.qpd_fields[field].component_index;
      for (unsigned int i = 0, i_field = 0; i_field < field_dofs_per_cell; /*increment at end of loop*/)
      {
        if (finite_element.system_to_component_index(i).first == field_component)
        {
          data.local_dof_indices[field][i_field] = scratch.local_dof_indices[i];
          ++i_field;
        }
        ++i;
      }

      // get the field values directly from QPDHandler.
      typename QPDHandler<dim>::active_cell_iterator qpd_cell(&triangulation,
                                                              cell->level(),
                                                              cell->index(),
                                                              &qpd_handler);

      for (unsigned int field = 0; field < n_fields; ++field)
        qpd_cell->get(scratch.qpd_fields[field].qpd_indicator, 
                      scratch.old_field_values[field]);

      if (parameters.use_ALE_method)
      {
        scratch.fe_values[introspection.extractors.displacement].get_function_values(
          solution, scratch.displacement_increments);
        scratch.fe_values[introspection.extractors.displacement].get_function_values(
          mesh_deformation_handler.mesh_displacement_increments, scratch.mesh_displacement_increments);
      }

      data.local_matrix = 0;
      for (unsigned int field = 0; field < data.local_rhs.size(); ++field)
        data.local_rhs[field] = 0;

      // trigger the invocation of the various functions that actually do
      // all of the assembling
      for (unsigned int i = 0; i < assemblers->qpd_system.size(); ++i)
        assemblers->qpd_system[i]->execute(scratch, data);

      // then also work on possible face terms.
      if (parameters.use_ALE_method)
      {
        // reset the matrices for interior face contributions
        for (auto &m : data.local_matrices_int_ext)
          m = 0;
        for (auto &m : data.local_matrices_ext_int)
          m = 0;
        for (auto &m : data.local_matrices_ext_ext)
          m = 0;

        // mark the arrays initialized to zero above as currently all unused.
        std::fill(data.assembled_matrices.begin(), data.assembled_matrices.end(), false);

        for (const unsigned int face_no : cell->face_indices())
        {
          scratch.face_number = face_no;
         
          (*scratch.fe_face_values).reinit(cell, face_no);
          (*scratch.fe_face_values)[introspection.extractors.displacement].get_function_values(
            solution, scratch.face_displacement_increments);
          (*scratch.fe_face_values)[introspection.extractors.displacement].get_function_values(
            mesh_deformation_handler.mesh_displacement_increments, scratch.face_mesh_displacement_increments);

          if (cell->at_boundary(face_no) && !cell->has_periodic_neighbor(face_no))
            for (unsigned int i = 0; i < assemblers->qpd_system_on_boundary_face.size(); ++i)
              assemblers->qpd_system_on_boundary_face[i]->execute(scratch, data);
          else
            for (unsigned int i = 0; i < assemblers->qpd_system_on_interior_face.size(); ++i)
              assemblers->qpd_system_on_interior_face[i]->execute(scratch, data);
        }
      }
    }
  }


  template <int dim>
  void
  Simulator<dim>::
  copy_local_to_global_qpd_system(const internal::Assembly::CopyData::QPDSystem<dim> &data)
  {
    // copy entries into the global matrix.
    current_constraints.distribute_local_to_global(data.local_matrix,
                                                   data.local_dof_indices[0],
                                                   system_matrix);

    // copy entries into the global RHS vector for each QPD field.
    for (unsigned int field = 0; field < data.local_rhs.size(); ++field)
      current_constraints.distribute_local_to_global(data.local_rhs[field],
                                                     data.local_dof_indices[field],
                                                     system_rhs);

    // copy DG contributions into the global matrix.
    if (parameters.use_ALE_method)
    {
      for (unsigned int f = 0; f < data.assembled_matrices.size(); ++f)
        if (data.assembled_matrices[f])
        {
          for (unsigned int i = 0; i < data.local_dof_indices[0].size(); ++i)
            for (unsigned int j = 0; j < data.neighbor_dof_indices[f].size(); ++j)
            {
              system_matrix.add(data.local_dof_indices[0][i],
                                data.neighbor_dof_indices[f][j],
                                data.local_matrices_int_ext[f](i, j));
              system_matrix.add(data.neighbor_dof_indices[f][j],
                                data.local_dof_indices[0][i],
                                data.local_matrices_ext_int[f](j, i));
            }

          for (unsigned int i = 0; i < data.neighbor_dof_indices[f].size(); ++i)
            for (unsigned int j = 0; j < data.neighbor_dof_indices[f].size(); ++j)
              system_matrix.add(data.neighbor_dof_indices[f][i],
                                data.neighbor_dof_indices[f][j],
                                data.local_matrices_ext_ext[f](i, j));
        }
    }
  }


  template <int dim>
  void Simulator<dim>::assemble_qpd_system(const std::vector<QPDField> &qpd_fields)
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

    // All QPD fields in the same category use the same block (the block of the
    // first field in the category) of system matrix for solving.
    const unsigned int block0_idx = qpd_fields[0].block_index;
    system_matrix.block(block0_idx, block0_idx) = 0;

    for (unsigned int field = 0; field < qpd_fields.size(); ++field)
      system_rhs.block(qpd_fields[field].block_index) = 0;

    using CellFilter = FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>;

    const UpdateFlags update_flags = update_values |
                                     update_JxW_values |
                                     (parameters.use_ALE_method 
                                      ?
                                      update_gradients
                                      :
                                      update_default);

    const UpdateFlags face_update_flags = (parameters.use_ALE_method 
                                           ?
                                           update_values |
                                           update_quadrature_points |
                                           update_normal_vectors |
                                           update_JxW_values
                                           :
                                           update_default);

    const unsigned int field_dofs_per_cell = finite_element.base_element(qpd_fields[0].base_index).dofs_per_cell;

    auto worker = [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
                      internal::Assembly::Scratch::QPDSystem<dim> &scratch,
                      internal::Assembly::CopyData::QPDSystem<dim> &data)
    {
      this->local_assemble_qpd_system(cell, scratch, data);
    };

    auto copier = [&](const internal::Assembly::CopyData::QPDSystem<dim> &data)
    {
      this->copy_local_to_global_qpd_system(data);
    };

    WorkStream::
    run(CellFilter(IteratorFilters::LocallyOwnedCell(),
                   dof_handler.begin_active()),
        CellFilter(IteratorFilters::LocallyOwnedCell(),
                   dof_handler.end()),
        worker,
        copier,
        internal::Assembly::Scratch::
        QPDSystem<dim>(finite_element,
                       *mapping,
                       quadrature_formula,
                       // Only generate a valid face quadrature if necessary
                       (parameters.use_ALE_method ?
                        QGauss<dim-1>(parameters.n_gaussian_points) :
                        Quadrature<dim-1>()),
                       update_flags,
                       face_update_flags,
                       field_dofs_per_cell,
                       qpd_fields),
        internal::Assembly::CopyData::
        QPDSystem<dim>(field_dofs_per_cell,
                       parameters.use_ALE_method,
                       qpd_fields));

    system_matrix.block(block0_idx, block0_idx).compress(VectorOperation::add);
    for (unsigned int field = 0; field < qpd_fields.size(); ++field)
      system_rhs.block(qpd_fields[field].block_index).compress(VectorOperation::add);
  }
}


// explicit instantiations
namespace elaspect
{
#define INSTANTIATE(dim) \
  template void Simulator<dim>::set_assemblers(); \
  template void Simulator<dim>::local_assemble_mechanical_system(const typename DoFHandler<dim>::active_cell_iterator &, \
                                                            internal::Assembly::Scratch::MechanicalSystem<dim> &, \
                                                            internal::Assembly::CopyData::MechanicalSystem<dim> &); \
  template void Simulator<dim>::copy_local_to_global_mechanical_system(const internal::Assembly::CopyData::MechanicalSystem<dim> &); \
  template void Simulator<dim>::assemble_mechanical_system(); \
  template void Simulator<dim>::local_assemble_thermo_system(const typename DoFHandler<dim>::active_cell_iterator &, \
                                                              internal::Assembly::Scratch::ThermoSystem<dim> &, \
                                                              internal::Assembly::CopyData::ThermoSystem<dim> &); \
  template void Simulator<dim>::copy_local_to_global_thermo_system(const internal::Assembly::CopyData::ThermoSystem<dim> &); \
  template void Simulator<dim>::assemble_thermo_system(); \
  template void Simulator<dim>::local_assemble_qpd_system(const typename DoFHandler<dim>::active_cell_iterator &, \
                                                          internal::Assembly::Scratch::QPDSystem<dim> &, \
                                                          internal::Assembly::CopyData::QPDSystem<dim> &); \
  template void Simulator<dim>::copy_local_to_global_qpd_system(const internal::Assembly::CopyData::QPDSystem<dim> &); \
  template void Simulator<dim>::assemble_qpd_system(const std::vector<QPDField> &);

  ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
