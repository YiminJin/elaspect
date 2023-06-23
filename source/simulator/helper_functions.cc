#include <elaspect/simulator.h>

#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/physics/elasticity/kinematics.h>
#include <deal.II/physics/transformations.h>

namespace elaspect
{
  namespace internal
  {
    template <int dim>
    unsigned int
    get_qpd_indicator(const Introspection<dim> &introspection,
                      const std::string        &name,
                      const unsigned int        component)
    {
      if (name == "composition")
      {
        AssertIndexRange(component, introspection.n_compositional_fields);
        return introspection.qpd_indicators.composition + component;
      }

      if (name == "stress")
      {
        AssertIndexRange(component, introspection.component_indices.stress.size());
        return introspection.qpd_indicators.old_stress + component;
      }

      if (name == "plastic strain")
      {
        AssertDimension(component, 0);
        return introspection.qpd_indicators.current_plastic_strain;
      }

      Assert(false, ExcNotImplemented());
      return numbers::invalid_unsigned_int;
    }
  }


  template <int dim>
  Simulator<dim>::QPDField::
  QPDField(const Introspection<dim> &introspection,
           const std::string        &name,
           const unsigned int        component)
    :
    field_type(name == "composition" ? compositional_field : physical_field),
    component_index(introspection.variable(name).first_component_index + component),
    block_index(introspection.variable(name).block_index + component),
    base_index(introspection.variable(name).base_index),
    qpd_indicator(internal::get_qpd_indicator(introspection, name, component)),
    compositional_variable(name == "composition" ? component : numbers::invalid_unsigned_int),
    scalar_extractor(FEValuesExtractors::Scalar(component_index))
  {}


  template <int dim>
  bool Simulator<dim>::QPDField::is_compositional_field() const
  {
    return (field_type == compositional_field);
  }


  template <int dim>
  void Simulator<dim>::output_statistics ()
  {
    // Only write the statistics file from processor zero
    if (Utilities::MPI::this_mpi_process(mpi_communicator) != 0)
      return;

    // Formatting the table we're about to output and writing the
    // actual file may take some time, so do it on a separate
    // thread. We do this using a lambda function that takes 
    // a copy of the statistics object to make sure that whatever
    // we do to the 'real' statistics object at the time of 
    // writing data doesn't affect what we write.
    //
    // Before we can start working on a new thread, we need to
    // make sure that the previous thread is done or they'll
    // step on each other's feet.
    if (output_statistics_thread.joinable())
      output_statistics_thread.join(); 

    // Write data in the background through a lambda function.
// This happening in the background means that we have
    // to create a copy of the statistics table, since whatever is
    // running in the foreground may continue to add entries to the
    // statistics table at the same time.
    output_statistics_thread = std::thread(
      [statistics_copy_ptr = std::make_unique<TableHandler>(statistics), this]()
{
        // First write everything into a string in memory
        std::ostringstream stream;
        statistics_copy_ptr->write_text (stream,
                                         TableHandler::table_with_separate_column_description);
        stream.flush();

        const std::string statistics_contents = stream.str();

        // Next find out whether we need to write everything into
        // the statistics file, or whether it is enough to just write
        // the last few bytes that were added since we wrote to that
        // file again. The way we do that is by checking whether the
        // first few bytes of the string we just created match what we
// had previously written. One might think that they always should,
        // but the statistics object automatically sizes the column widths
        // of its output to match what is being written, and so if a later
        // entry requires more width, then even the first columns are
        // changed -- in that case, we will have to write everything,
        // not just append one line.
        const bool write_everything
          = ( // We may have never written anything. More precisely, this
              // case happens if the statistics_last_write_size is at the
              // value initialized by the Simulator::Simulator()
              // constructor, and this can happen in two situations:
              // (i) At the end of the first time step; and (ii) upon restart
              // since the variable we query here is not serialized. It is clear
              // that in both situations, we want to write the
              // entire contents of the statistics object. For the second
              // case, this is also appropriate since someone may have
              // previously restarted from a checkpoint, run a couple of
              // time steps that have added to the statistics file, but then
              // aborted the run again; a later restart from the same
              // checkpoint then requires overwriting the statistics file
              // on disk with what we have when this function is called for
              // the first time after the restart. The same situation
              // happens if the simulation kept running for some time after
              // a checkpoint, but is resumed from that checkpoint (i.e.,
              // at an earlier time step than when the statistics file was
              // written to last). In these situations, we effectively want
              // to "truncate" the file to the state stored in the checkpoint,
              // and we do that by just overwriting the entire file.
              (statistics_last_write_size == 0)
              ||
              // Or the size of the statistics file may have
              // shrunk mysteriously -- this shouldn't happen
              // but if it did we'd get into trouble with the
              // .substr() call in the next check.
              (statistics_last_write_size > statistics_contents.size())
              ||
              // Or the hash of what we wrote last time doesn't match
              // the hash of the first part of what we want to write
              (statistics_last_hash
               !=
               std::hash<std::string>()(statistics_contents.substr(0, statistics_last_write_size))) );

        const std::string stat_file_name = parameters.output_directory + "statistics";

        if (write_everything)
        {
          // Write what we have into a tmp file, then move that into
          // place
const std::string tmp_file_name = stat_file_name + ".tmp";
          {
            std::ofstream tmp_file (tmp_file_name);
            tmp_file << statistics_contents;
          }
          std::rename(tmp_file_name.c_str(), stat_file_name.c_str());
        }
        else
        {
          // If we don't have to write everything, then the first part of what
          // we want to write matches what's already on disk. In that case,
          // we just have to append what's new.
          std::ofstream stat_file (stat_file_name, std::ios::app);
          stat_file << statistics_contents.substr(statistics_last_write_size, std::string::npos);
        }

        // Now update the size and hash of what we just wrote so that
        // we can compare against it next time we get here. Note that we do
        // not need to guard access to these variables with a mutex because
        // this is the only function that touches the variables, and
        // this function runs only once at a time (on a different
        // thread, but it's not started a second time while the previous
        // run hasn't finished).
        statistics_last_write_size = statistics_contents.size();
        statistics_last_hash       = std::hash<std::string>()(statistics_contents);
      });
  }

  template <int dim>
  bool 
  Simulator<dim>::
  maybe_do_initial_refinement (const unsigned int max_refinement_level)
  {
    if (pre_refinement_step < parameters.initial_adaptive_refinement)
    {
      if (parameters.timing_output_frequency == 0)
        computing_timer.print_summary();

      output_statistics();

      if (parameters.run_postprocessors_on_initial_refinement &&
          (!parameters.run_postprocessors_on_nonlinear_iterations))
        postprocess();

      refine_mesh (max_refinement_level);
      ++pre_refinement_step;
      return true;
    }
    else
    {
      // invalidate the value of pre_refinement_step since it will no longer be used from here on
      pre_refinement_step = std::numeric_limits<unsigned int>::max();
      return false;
    }
  }


  template <int dim>
  void 
  Simulator<dim>::maybe_refine_mesh (const unsigned int max_refinement_level)
  {
    if ((timestep_number > 0
         &&
         parameters.adaptive_refinement_interval > 0
         &&
         timestep_number % parameters.adaptive_refinement_interval == 0)
        ||
        (timestep_number == 0 && parameters.adaptive_refinement_interval == 1))
    {
      refine_mesh (max_refinement_level);
    }
  }


  template <int dim>
  void Simulator<dim>::maybe_write_timing_output () const
  {
    bool write_timing_output = false;
    if (parameters.timing_output_frequency <= 1)
      write_timing_output = true;
    else if ((timestep_number > 0) &&
             (timestep_number % parameters.timing_output_frequency == 0))
      write_timing_output = true;

    // if requested output a summary of the current timing information
    if (write_timing_output)
      computing_timer.print_summary();
  }


  template <int dim>
  void 
  Simulator<dim>::
  interpolate_onto_displacement_system (const TensorFunction<1,dim> &func,
                                        TrilinosWrappers::MPI::Vector &vec)
  {
    AffineConstraints<double> hanging_constraints(introspection.index_sets.system_relevant_set);
    DoFTools::make_hanging_node_constraints(dof_handler, hanging_constraints);
    hanging_constraints.close();

    const std::vector<Point<dim> > mesh_support_points = finite_element.base_element(introspection.base_elements.displacement).get_unit_support_points();
    FEValues<dim> mesh_points (*mapping, finite_element, mesh_support_points, update_quadrature_points);
    std::vector<types::global_dof_index> cell_dof_indices (finite_element.dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
      {
        mesh_points.reinit(cell);
        cell->get_dof_indices (cell_dof_indices);
        for (unsigned int j=0; j<finite_element.base_element(introspection.base_elements.displacement).dofs_per_cell; ++j)
          for (unsigned int dir=0; dir<dim; ++dir)
          {
            unsigned int support_point_index
              = finite_element.component_to_system_index(introspection.component_indices.displacement[dir],
                                                                                 /*dof index within component=*/ j);
            Assert(introspection.block_indices.displacement == 0, ExcNotImplemented());
            vec[cell_dof_indices[support_point_index]] = func.value(mesh_points.quadrature_point(j))[dir];
          }
      }

    vec.compress(VectorOperation::insert);
    hanging_constraints.distribute(vec);
  }


  template <int dim>
  void
  Simulator<dim>::apply_return_mapping()
  {
    AssertThrow(parameters.constitutive_relation & ConstitutiveRelation::plasticity,
                ExcInternalError());

    FEValues<dim> fe_values(*mapping, 
                            finite_element,
                            quadrature_formula,
                            update_values | 
                            update_gradients | 
                            update_quadrature_points);

    const unsigned int n_q_points = fe_values.n_quadrature_points;

    const MaterialModel::FieldDependences::Dependence field_dependences
      = material_handler.get_field_dependences_for_return_mapping();
    const MaterialModel::MaterialProperties::Property requested_properties
      = material_handler.get_material_properties_for_return_mapping();
    MaterialModel::MaterialModelInputs<dim> material_inputs(n_q_points,
                                                            introspection.n_compositional_fields,
                                                            field_dependences,
                                                            requested_properties);

    std::vector<SymmetricTensor<2,dim>> stresses(n_q_points);
    std::vector<double> plastic_strains(n_q_points);
    std::vector<double> flags(n_q_points);

    const unsigned int stress_indicator         = introspection.qpd_indicators.current_stress;
    const unsigned int plastic_strain_indicator = introspection.qpd_indicators.current_plastic_strain;
    const unsigned int flag_indicator           = introspection.qpd_indicators.return_mapping_flag;

    typename DoFHandler<dim>::active_cell_iterator
    dof_cell = dof_handler.begin_active(),
    dof_endc = dof_handler.end();
    typename QPDHandler<dim>::active_cell_iterator
    qpd_cell = qpd_handler.begin_active();
    for (; dof_cell != dof_endc; ++dof_cell, ++qpd_cell)
      if (dof_cell->is_locally_owned())
      {
        fe_values.reinit(dof_cell);
        material_inputs.reinit(fe_values, qpd_handler, introspection, solution);

        material_handler.apply_return_mapping(material_inputs, 
                                              stresses,
                                              plastic_strains,
                                              flags);

        for (unsigned int q = 0; q < n_q_points; ++q)
          qpd_cell->set(q, stress_indicator, stresses[q]);

        qpd_cell->set(plastic_strain_indicator, plastic_strains);
        qpd_cell->set(flag_indicator, flags);
      }
  }


  namespace
  {
    Tensor<2,2> get_rotation_matrix (const Tensor<2,2> &grad_u)
    {
      const double curl = grad_u[1][0] - grad_u[0][1];
      const double angle = std::atan(curl);
      return dealii::Physics::Transformations::Rotations::rotation_matrix_2d(-angle);
    }


    Tensor<2,3> get_rotation_matrix (const Tensor<2,3> &grad_u)
    {
      const Point<3> curl (grad_u[2][1] - grad_u[1][2], 
                           grad_u[0][2] - grad_u[2][0],
                           grad_u[1][0] - grad_u[0][1]);

      const double tan_angle = std::sqrt(curl * curl);
      const double angle     = std::atan(tan_angle);

      if (std::abs(angle) < 1e-9)
      {
        static const double rotation[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        static const Tensor<2,3> rot(rotation);
        return rot;
      }

      const Point<3> axis = curl / tan_angle;
      return Physics::Transformations::Rotations::rotation_matrix_3d(axis, -angle);
    }

  }


  template <int dim>
  void Simulator<dim>::update_quadrature_point_data()
  {
    FEValues<dim> fe_values(*mapping, finite_element,
                            quadrature_formula,
                            update_values | 
                            update_gradients | 
                            update_quadrature_points);

    MaterialModel::FieldDependences::Dependence dependences = MaterialModel::FieldDependences::displacement_gradient;
    MaterialModel::MaterialProperties::Property properties = MaterialModel::MaterialProperties::reaction_terms;

    if (!(parameters.constitutive_relation & ConstitutiveRelation::plasticity))
    {
      properties |= MaterialModel::MaterialProperties::tangent_modulus;

      if (parameters.constitutive_relation & ConstitutiveRelation::viscosity)
        properties |= MaterialModel::MaterialProperties::viscosity;

      if (parameters.constitutive_relation & ConstitutiveRelation::thermal_expansion)
      {
        dependences |= MaterialModel::FieldDependences::temperature;
        properties |= MaterialModel::MaterialProperties::thermal_expansivity;
      }
    }
    
    dependences |= material_handler.get_field_dependences_for_evaluation(properties);

    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int n_comp = introspection.n_compositional_fields;
    MaterialModel::MaterialModelInputs<dim> material_inputs(n_q_points, n_comp, dependences, properties);
    MaterialModel::MaterialModelOutputs<dim> material_outputs(n_q_points, n_comp);

    std::vector<double> qpd_component_values(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator
    dof_cell = dof_handler.begin_active(),
    dof_endc = dof_handler.end();
    typename QPDHandler<dim>::active_cell_iterator
    qpd_cell = qpd_handler.begin_active();
    for (; dof_cell != dof_endc; ++dof_cell, ++qpd_cell)
      if (dof_cell->is_locally_owned())
      {
        fe_values.reinit(dof_cell);
        material_inputs.reinit(fe_values, qpd_handler, introspection, solution);
        material_handler.evaluate(material_inputs, material_outputs);

        // deal with compositional fields:
        for (unsigned int i = 0; i < n_comp; ++i)
        {
          qpd_cell->get(introspection.qpd_indicators.composition + i, qpd_component_values);
          for (unsigned int q = 0; q < n_q_points; ++q)
            qpd_component_values[q] += material_outputs.reaction_terms[q][i];
          qpd_cell->set(introspection.qpd_indicators.composition + i, qpd_component_values);
        }

        // deal with stress field:
        for (unsigned int q = 0; q < n_q_points; ++q)
        {
          SymmetricTensor<2,dim> stress;
          if (parameters.constitutive_relation & ConstitutiveRelation::plasticity)
            // If plasticity is included in the constitutive relationship, then stress has 
            // already been updated in the return-mapping procedure.
            stress = qpd_cell->get_symmetric_tensor(q, introspection.qpd_indicators.current_stress);
          else
          {
            // Otherwise, the stress should be updated here through the incremental equation
            // \sigma = \sigma_{old} + D : \Delta\varepsilon.
            stress = qpd_cell->get_symmetric_tensor(q, introspection.qpd_indicators.old_stress);
            stress += material_outputs.tangent_moduli[q] * 
                      symmetrize(material_inputs.incremental_displacement_gradient[q]);

            // eliminate viscous relaxation
            if (parameters.constitutive_relation & ConstitutiveRelation::viscosity)
            {
              double Gve = MaterialUtilities::get_shear_modulus(material_outputs.tangent_moduli[q]);
              double eta = material_outputs.viscosities[q];

              double pressure = trace(stress) / dim;
              SymmetricTensor<2,dim> deviatoric_stress = deviator(stress) * ((eta - Gve * time_step) / eta);

              stress = deviatoric_stress + pressure * unit_symmetric_tensor<dim>();
            }

            // eliminate thermal expansion
            if (parameters.constitutive_relation & ConstitutiveRelation::thermal_expansion)
            {
              const double T_old = material_inputs.qpd_cell->get_scalar(
                q, introspection.qpd_indicators.old_temperature);

              const double alpha_dT = 
                (timestep_number > 0 ?
                 material_outputs.thermal_expansivities[q] * (material_inputs.temperature[q] - T_old) :
                 0.0);

              stress -= material_outputs.tangent_moduli[q] *
                        (alpha_dT * unit_symmetric_tensor<dim>());
            }
          }

          // rotate the stress tensor
          const Tensor<2,dim> rotation = 
            get_rotation_matrix(material_inputs.incremental_displacement_gradient[q]);
          const SymmetricTensor<2,dim> rotated_stress = 
            symmetrize(transpose(rotation) * static_cast<Tensor<2,dim>>(stress) * rotation);

          qpd_cell->set(q, introspection.qpd_indicators.old_stress, rotated_stress);
          if (parameters.constitutive_relation & ConstitutiveRelation::plasticity)
            qpd_cell->set(q, introspection.qpd_indicators.current_stress, rotated_stress);
        }

        // If plasticity is included in the constitutive relationship, then 
        // update the values of old plastic strain and set return-mapping flags
        // to zero.
        if (parameters.constitutive_relation & ConstitutiveRelation::plasticity)
        {
          qpd_cell->get(introspection.qpd_indicators.current_plastic_strain, qpd_component_values);
          qpd_cell->set(introspection.qpd_indicators.old_plastic_strain, qpd_component_values);

          std::fill(qpd_component_values.begin(), qpd_component_values.end(), 0.0);
          qpd_cell->set(introspection.qpd_indicators.return_mapping_flag, qpd_component_values);
        }
      }
  }


  template <int dim>
  void
  Simulator<dim>::advect_quadrature_point_data()
  {
    // If ALE method is applied, then assemble and solve a convection equation
    // for each QPD component; otherwise, just project the QPD components onto
    // grid nodes.
    const std::string section_name = 
      (parameters.use_ALE_method ? "QPD advection" : "QPD smoothing");
    TimerOutput::Scope timer(computing_timer, section_name);

    UpdateFlags update_flags(update_values | update_JxW_values);
    if (parameters.use_ALE_method)
      update_flags |= update_gradients;

    FEValues<dim> fe_values(*mapping, finite_element, quadrature_formula, update_flags);

    const unsigned int n_q_points = fe_values.n_quadrature_points;
    const unsigned int dofs_per_cell = fe_values.dofs_per_cell;

    // get the component index, block index and data indicator for each QPD component
    std::vector<unsigned int> qpd_components, qpd_blocks, qpd_indicators;
    for (unsigned int i = 0; i < introspection.n_compositional_fields; ++i)
    {
      qpd_components.push_back(introspection.component_indices.compositional_fields[i]);
      qpd_blocks.push_back(introspection.block_indices.compositional_fields[i]);
      qpd_indicators.push_back(introspection.qpd_indicators.composition + i);
    }
    for (unsigned int i = 0; i < SymmetricTensor<2,dim>::n_independent_components; ++i)
    {
      qpd_components.push_back(introspection.component_indices.stress[i]);
      qpd_blocks.push_back(introspection.block_indices.stress[i]);
      qpd_indicators.push_back(introspection.qpd_indicators.old_stress + i);
    }
    if (parameters.constitutive_relation & ConstitutiveRelation::plasticity)
    {
      const FEVariable<dim> &variable = introspection.variable("plastic strain");
      qpd_components.push_back(variable.first_component_index);
      qpd_blocks.push_back(variable.block_index);
      qpd_indicators.push_back(introspection.qpd_indicators.current_plastic_strain);
    }

    // notice that compositional fields and physical fields (stress, plastic strain, etc.) 
    // may have different polynomial degree, so we need to assemble their matrices
    // separatedly. The first component of compositional/physical fields is used to 
    // calculate the shape functions
    const unsigned int first_phys_component = introspection.n_compositional_fields;

    const FEValuesExtractors::Scalar comp_extractor(qpd_components[0]);
    const FEValuesExtractors::Scalar phys_extractor(qpd_components[first_phys_component]);

    const unsigned int comp_dofs_per_cell = finite_element.base_element(introspection.base_elements.compositional_fields).dofs_per_cell;
    const unsigned int phys_dofs_per_cell = finite_element.base_element(introspection.base_elements.stress).dofs_per_cell;

    FullMatrix<double> cell_comp_matrix(comp_dofs_per_cell, comp_dofs_per_cell);
    FullMatrix<double> cell_phys_matrix(phys_dofs_per_cell, phys_dofs_per_cell);

    std::vector<double> phi_comp(comp_dofs_per_cell);
    std::vector<double> phi_phys(phys_dofs_per_cell);

    // convection term: u\dot\nabla\phi
    std::vector<double> u_dot_grad_phi_comp(comp_dofs_per_cell);
    std::vector<double> u_dot_grad_phi_phys(comp_dofs_per_cell);

    std::vector<types::global_dof_index> cell_dof_indices(dofs_per_cell);

    std::vector<Vector<double>> cell_rhs;
    std::vector<std::vector<types::global_dof_index>> qpd_dof_indices;
    for (unsigned int i = 0; i < qpd_components.size(); ++i)
    {
      if (i < introspection.n_compositional_fields)
      {
        cell_rhs.push_back(Vector<double>(comp_dofs_per_cell));
        qpd_dof_indices.push_back(std::vector<types::global_dof_index>(comp_dofs_per_cell));
      }
      else
      {
        cell_rhs.push_back(Vector<double>(phys_dofs_per_cell));
        qpd_dof_indices.push_back(std::vector<types::global_dof_index>(phys_dofs_per_cell));
      }
    }

    std::vector<std::vector<double>> qpd_component_values(
      qpd_components.size(), std::vector<double>(n_q_points));

    // we use the first block of compositional/physical fields to store the matrix
    TrilinosWrappers::SparseMatrix &comp_matrix = system_matrix.block(qpd_blocks[0], qpd_blocks[0]);
    TrilinosWrappers::SparseMatrix &phys_matrix = system_matrix.block(qpd_blocks[first_phys_component], 
                                                                      qpd_blocks[first_phys_component]);

    // data for calculating the relative displacement increment
    std::vector<Tensor<1,dim>> displacement_increments(n_q_points);
    std::vector<Tensor<1,dim>> mesh_displacement_increments(n_q_points);
    std::vector<Tensor<1,dim>> relative_displacement_increments(n_q_points);

    const bool has_compositional_fields = (introspection.n_compositional_fields > 0);

    typename DoFHandler<dim>::active_cell_iterator
    dof_cell = dof_handler.begin_active(),
    dof_endc = dof_handler.end();
    typename QPDHandler<dim>::active_cell_iterator
    qpd_cell = qpd_handler.begin_active();
    for (; dof_cell != dof_endc; ++dof_cell, ++qpd_cell)
      if (dof_cell->is_locally_owned())
      {
        fe_values.reinit(dof_cell);
        dof_cell->get_dof_indices(cell_dof_indices);

        // sort out the dof indices for each QPD component
        for (unsigned int comp = 0; comp < qpd_components.size(); ++comp)
        {
          qpd_cell->get(qpd_indicators[comp], qpd_component_values[comp]);

          for (unsigned int i = 0, i_qpd = 0; i_qpd < qpd_dof_indices[comp].size(); /*increment at end of loop*/)
          {
            if (finite_element.system_to_component_index(i).first == qpd_components[comp])
            {
              qpd_dof_indices[comp][i_qpd] = cell_dof_indices[i];
              ++i_qpd;
            }
            ++i;
          }
        }

        if (parameters.use_ALE_method)
        {
          fe_values[introspection.extractors.displacement].get_function_values(solution, displacement_increments);
          fe_values[introspection.extractors.displacement].get_function_values(mesh_deformation_handler.mesh_displacement_increments, 
                                                                               mesh_displacement_increments);
        } 

        cell_comp_matrix = 0;
        cell_phys_matrix = 0;
        for (unsigned int comp = 0; comp < qpd_components.size(); ++comp)
          cell_rhs[comp] = 0;

        // calculate the SUPG parameter k
        double k_comp = 0, k_phys = 0;
        if (parameters.use_ALE_method)
        {
          double u_max = 0;
          for (unsigned int q = 0; q < n_q_points; ++q)
          {
            relative_displacement_increments[q] = 
              displacement_increments[q] - mesh_displacement_increments[q];
            u_max = std::max(relative_displacement_increments[q].norm(), u_max);
          }

          double k = 0;
          if(timestep_number > 0 && u_max > 1e-50)
            k = dof_cell->diameter() / (2. * u_max);

          if (has_compositional_fields)
            k_comp = k / finite_element.base_element(introspection.base_elements.compositional_fields).degree;
          k_phys = k / finite_element.base_element(introspection.base_elements.stress).degree;
        }

        for (unsigned int q = 0; q < n_q_points; ++q)
        {
          // calculate the shape functions of the first component of 
          // compositional/physical fields
          if (has_compositional_fields)
          {
            for (unsigned int i = 0, i_comp = 0; i_comp < comp_dofs_per_cell; /*increment at end of loop*/)
            {
              if (finite_element.system_to_component_index(i).first == qpd_components[0])
              {
                phi_comp[i_comp] = fe_values[comp_extractor].value(i, q);
                if (parameters.use_ALE_method)
                  u_dot_grad_phi_comp[i_comp] = relative_displacement_increments[q] * 
                                                fe_values[comp_extractor].gradient(i, q);

                ++i_comp;
              }
              ++i;
            }
          }

          for (unsigned int i = 0, i_phys = 0; i_phys < phys_dofs_per_cell; /*increment at end of loop*/)
          {
            if (finite_element.system_to_component_index(i).first == qpd_components[first_phys_component])
            {
              phi_phys[i_phys] = fe_values[phys_extractor].value(i, q);
              if (parameters.use_ALE_method)
                u_dot_grad_phi_phys[i_phys] = relative_displacement_increments[q] * 
                                              fe_values[phys_extractor].gradient(i, q);

              ++i_phys;
            }
            ++i;
          }

          const double JxW = fe_values.JxW(q);

          // assemble the cell matrix for compositional/physical field
          if (parameters.use_ALE_method)
          {
            if (has_compositional_fields)
            {
              for (unsigned int i = 0; i < comp_dofs_per_cell; ++i)
                for (unsigned int j = 0; j < comp_dofs_per_cell; ++j)
                  cell_comp_matrix(i, j) += (phi_comp[i] + k_comp * u_dot_grad_phi_comp[i])
                                            * (phi_comp[j] + u_dot_grad_phi_comp[j])
                                            * JxW;
            }
            
            for (unsigned int i = 0; i < phys_dofs_per_cell; ++i)
              for (unsigned int j = 0; j < phys_dofs_per_cell; ++j)
                cell_phys_matrix(i, j) += (phi_phys[i] + k_phys * u_dot_grad_phi_phys[i])
                                          * (phi_phys[j] + u_dot_grad_phi_phys[j])
                                          * JxW;
          }
          else
          {
            if (has_compositional_fields)
            {
              for (unsigned int i = 0; i < comp_dofs_per_cell; ++i)
                for (unsigned int j = 0; j < comp_dofs_per_cell; ++j)
                  cell_comp_matrix(i, j) += (phi_comp[i] * phi_comp[j]) * JxW;
            }

            for (unsigned int i = 0; i < phys_dofs_per_cell; ++i)
              for (unsigned int j = 0; j < phys_dofs_per_cell; ++j)
                cell_phys_matrix(i, j) += (phi_phys[i] * phi_phys[j]) * JxW;
          }

          // assemble the rhs of each component of compositional/physical field
          for (unsigned int comp = 0; comp < qpd_components.size(); ++comp)
          {
            if (parameters.use_ALE_method)
            {
              if (comp < introspection.n_compositional_fields)
              {
                for (unsigned int i = 0; i < comp_dofs_per_cell; ++i)
                  cell_rhs[comp](i) += (phi_comp[i] + k_comp * u_dot_grad_phi_comp[i])
                                       * qpd_component_values[comp][q] * JxW;
              }
              else
              {
                for (unsigned int i = 0; i < phys_dofs_per_cell; ++i)
                  cell_rhs[comp](i) += (phi_phys[i] + k_phys * u_dot_grad_phi_phys[i])
                                       * qpd_component_values[comp][q] * JxW;
              }
            }
            else
            {
              if (comp < introspection.n_compositional_fields)
              {
                for (unsigned int i = 0; i < comp_dofs_per_cell; ++i)
                  cell_rhs[comp](i) += (phi_comp[i] * qpd_component_values[comp][q]) * JxW;
              }
              else
              {
                for (unsigned int i = 0; i < phys_dofs_per_cell; ++i)
                  cell_rhs[comp](i) += (phi_phys[i] * qpd_component_values[comp][q]) * JxW;
              }
            }
          }
        }

        if (has_compositional_fields)
          current_constraints.distribute_local_to_global(cell_comp_matrix,
                                                         qpd_dof_indices[0],
                                                         system_matrix);

        current_constraints.distribute_local_to_global(cell_phys_matrix,
                                                       qpd_dof_indices[first_phys_component],
                                                       system_matrix);

        for (unsigned int comp = 0; comp < qpd_components.size(); ++comp)
          current_constraints.distribute_local_to_global(cell_rhs[comp],
                                                         qpd_dof_indices[comp],
                                                         system_rhs);
      }

    if (has_compositional_fields)
      comp_matrix.compress(VectorOperation::add);
    phys_matrix.compress(VectorOperation::add);
    for (unsigned int comp = 0; comp < qpd_components.size(); ++comp)
      system_rhs.block(qpd_blocks[comp]).compress(VectorOperation::add);

    // solve for each QPD component
    TrilinosWrappers::PreconditionILU comp_preconditioner, phys_preconditioner;
    if (has_compositional_fields)
      comp_preconditioner.initialize(comp_matrix);
    phys_preconditioner.initialize(phys_matrix);

    TrilinosWrappers::MPI::BlockVector distributed_solution(
      introspection.index_sets.system_partitioning, mpi_communicator);
    for (unsigned int comp = 0; comp < qpd_components.size(); ++comp)
      distributed_solution.block(qpd_blocks[comp]) = solution.block(qpd_blocks[comp]);
    current_constraints.set_zero(distributed_solution);

    for (unsigned int comp = 0; comp < qpd_components.size(); ++comp)
    {
      const double tolerance = 1e-8 * system_rhs.block(qpd_blocks[comp]).l2_norm();
      if (tolerance < 1e-50)
      {
        distributed_solution.block(qpd_blocks[comp]) = 0;
        continue;
      }

      SolverControl solver_control(1000, tolerance);
      std::unique_ptr<TrilinosWrappers::SolverBase> solver;
      if (parameters.use_ALE_method)
        solver = std::make_unique<TrilinosWrappers::SolverGMRES>(solver_control);
      else
        solver = std::make_unique<TrilinosWrappers::SolverCG>(solver_control);

      if (comp < introspection.n_compositional_fields)
        solver->solve(comp_matrix,
                      distributed_solution.block(qpd_blocks[comp]),
                      system_rhs.block(qpd_blocks[comp]),
                      comp_preconditioner);
      else
        solver->solve(phys_matrix,
                      distributed_solution.block(qpd_blocks[comp]),
                      system_rhs.block(qpd_blocks[comp]),
                      phys_preconditioner);
    }

    current_constraints.distribute(distributed_solution);
    for (unsigned int comp = 0; comp < qpd_components.size(); ++comp)
      solution.block(qpd_blocks[comp]) = distributed_solution.block(qpd_blocks[comp]);

    // interpolate the QPD components from grid nodes to quadrature points
    refresh_quadrature_point_data();
  }


  template <int dim>
  void Simulator<dim>::refresh_quadrature_point_data()
  {
    // TODO: determine whether to interpolate the QPD components in different
    // situations. For instance, all compositional and physical quantities 
    // must be refreshed after mesh refinement; if the mesh stays unchanged 
    // and the model is in Lagrangian description, then compositional fields
    // and plastic strain can be left as they were so that no artificial 
    // diffusion would be introduced to these fields.
    FEValues<dim> fe_values(*mapping, finite_element, quadrature_formula, update_values);

    const unsigned int n_q_points = fe_values.n_quadrature_points;
    std::vector<double> qpd_component_values(n_q_points);
    std::vector<SymmetricTensor<2,dim>> stress_values(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator
    dof_cell = dof_handler.begin_active(),
    dof_endc = dof_handler.end();
    typename QPDHandler<dim>::active_cell_iterator
    qpd_cell = qpd_handler.begin_active();
    for (; dof_cell != dof_endc; ++dof_cell, ++qpd_cell)
      if (dof_cell->is_locally_owned())
      {
        fe_values.reinit(dof_cell);
        
        // interpolate compositional fields
        for (unsigned int i = 0; i < introspection.n_compositional_fields; ++i)
        {
          fe_values[introspection.extractors.compositional_fields[i]].get_function_values(
            solution, qpd_component_values);
          qpd_cell->set(introspection.qpd_indicators.composition + i, qpd_component_values);
        }

        // interpolate stress values
        fe_values[introspection.extractors.stress].get_function_values(solution, stress_values);
        for (unsigned int q = 0; q < n_q_points; ++q)
        {
          qpd_cell->set(q, introspection.qpd_indicators.old_stress, stress_values[q]);
          if (parameters.constitutive_relation & ConstitutiveRelation::plasticity)
            qpd_cell->set(q, introspection.qpd_indicators.current_stress, stress_values[q]);
        }

        // interpolate plastic strain values
        if (parameters.constitutive_relation & ConstitutiveRelation::plasticity)
        {
          const FEValuesExtractors::Scalar extractor(introspection.variable("plastic strain").first_component_index);
          fe_values[extractor].get_function_values(solution, qpd_component_values);
          qpd_cell->set(introspection.qpd_indicators.old_plastic_strain, qpd_component_values);
          qpd_cell->set(introspection.qpd_indicators.current_plastic_strain, qpd_component_values);
        }
      }
  }
}


// explicit instantiations
namespace elaspect
{
#define INSTANTIATE(dim) \
  template struct Simulator<dim>::QPDField; \
  template void Simulator<dim>::output_statistics(); \
  template bool Simulator<dim>::maybe_do_initial_refinement (const unsigned int max_refinement_level); \
  template void Simulator<dim>::maybe_refine_mesh (const unsigned int max_refinement_level); \
  template void Simulator<dim>::maybe_write_timing_output () const; \
  template void Simulator<dim>::interpolate_onto_displacement_system(const TensorFunction<1,dim> &func, \
                                                                     TrilinosWrappers::MPI::Vector &vec); \
  template void Simulator<dim>::apply_return_mapping(); \
  template void Simulator<dim>::update_quadrature_point_data(); \
  template void Simulator<dim>::advect_quadrature_point_data(); \
  template void Simulator<dim>::refresh_quadrature_point_data();

  ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
