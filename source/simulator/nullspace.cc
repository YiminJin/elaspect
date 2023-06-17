#include <elaspect/simulator.h>
#include <elaspect/geometry_model/sphere.h>

namespace elaspect
{
  namespace internal
  {
    /**
     * A class we use when setting up the data structures for nullspace removal
     * of the rotations in spherical or annular shells.
     */
    template <int dim>
    class Rotation : public TensorFunction<1,dim>
    {
      private:
        const Tensor<1,dim> axis;

      public:
        // Constructor for TensorFunction that takes cartesian direction (1,2, or 3)
        // and creates a solid body rotation around that axis.
        Rotation(const unsigned int a)
          :
          axis(Tensor<1,dim>(Point<dim>::unit_vector(a)))
        {}

        // Constructor for TensorFunction that takes an axis
        // and creates a solid body rotation around that axis.
        Rotation(const Tensor<1,dim> &rotation_axis)
          :
          axis(rotation_axis)
        {}

        Tensor<1,dim> value (const Point<dim> &p) const override
        {
          if ( dim == 2)
            return cross_product_2d(p);
          else
            return cross_product_3d(axis, p);
        }
    };


    /**
     * A class we use when setting up the data structures for nullspace removal
     * of the translations in box-like geometries.
     */
    template <int dim>
    class Translation : public TensorFunction<1,dim>
    {
      private:
        const Tensor<1,dim> translation;

      public:
        // Constructor for TensorFunction that takes a Cartesian direction (1,2, or 3)
        // and creates a translation along that axis
        Translation(const unsigned int d)
          :
          translation( Point<dim>::unit_vector(d) )
        {}

        // Constructor for TensorFunction that takes a vector
        // and creates a translation along that vector
        Translation(const Tensor<1,dim> &t)
          :
          translation(t)
        {}

        Tensor<1,dim> value(const Point<dim> &) const override
        {
          return translation;
        }
    };
  }


  template <int dim>
  void
  Simulator<dim>::
  setup_nullspace_constraints (AffineConstraints<double> &constraints)
  {
    if (!(parameters.nullspace_removal & (NullspaceRemoval::linear_momentum
                                          | NullspaceRemoval::net_translation)))
      return;

    // Note: We want to add a single Dirichlet zero constraint for each
    // translation direction. This is complicated by the fact that we need to
    // find a DoF that is not already constrained. In parallel the constraint
    // needs to be added on all processors where it is locally_relevant and
    // all processors need to agree on the index.

    // First find candidates for DoF indices to constrain for each displacement component.
    types::global_dof_index disp_idx[dim];
    {
      for (unsigned int d = 0; d < dim; ++d)
        disp_idx[d] = numbers::invalid_dof_index;

      unsigned int n_left_to_find = dim;

      std::vector<types::global_dof_index> local_dof_indices (finite_element.dofs_per_cell);
      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned())
        {
          cell->get_dof_indices (local_dof_indices);

          for (unsigned int i = 0; i < finite_element.dofs_per_cell; ++i)
          {
            const unsigned int component = finite_element.system_to_component_index(i).first;

            if (component < introspection.component_indices.displacement[0]
                || component > introspection.component_indices.displacement[dim-1])
              continue; // only look at displacement

            const unsigned int disp_component = component - introspection.component_indices.displacement[0];

            if (disp_idx[disp_component] != numbers::invalid_dof_index)
              continue; // already found one

            const types::global_dof_index idx = local_dof_indices[i];

            if (constraints.can_store_line(idx) && !constraints.is_constrained(idx))
            {
              disp_idx[disp_component] = idx;
              --n_left_to_find;
            }

            // are we done searching?
            if (n_left_to_find == 0)
              break; // exit inner loop
          }

          if (n_left_to_find == 0)
            break; // exit outer loop
        }
    }

    const unsigned int flags[] = {(NullspaceRemoval::linear_momentum_x
                                   |NullspaceRemoval::net_translation_x),
                                  (NullspaceRemoval::linear_momentum_y
                                   |NullspaceRemoval::net_translation_y),
                                  (NullspaceRemoval::linear_momentum_z
                                   |NullspaceRemoval::net_translation_z)
                                 };

    for (unsigned int d=0; d<dim; ++d)
      if (parameters.nullspace_removal & flags[d])
      {
        // Make a reduction to find the smallest index (processors that
        // found a larger candidate just happened to not be able to store
        // that index with the minimum value). Note that it is possible that
        // some processors might not be able to find a potential DoF, for
        // example because they don't own any DoFs. On those processors we
        // will use dof_handler.n_dofs() when building the minimum (larger
        // than any valid DoF index).
        const types::global_dof_index global_idx = dealii::Utilities::MPI::min(
                                                     (disp_idx[d] != numbers::invalid_dof_index)
                                                     ?
                                                     disp_idx[d]
                                                     :
                                                     dof_handler.n_dofs(),
                                                     mpi_communicator);

        Assert(global_idx < dof_handler.n_dofs(),
               ExcMessage("Error, couldn't find a displacement DoF to constrain."));

        // Finally set this DoF to zero (if the current MPI process
        // cares about it):
        if (constraints.can_store_line(global_idx))
        {
          Assert(!constraints.is_constrained((global_idx)),
                 ExcInternalError());
          constraints.add_line(global_idx);
        }
      }
  }


  template <int dim>
  void 
  Simulator<dim>::
  remove_nullspace(TrilinosWrappers::MPI::BlockVector &relevant_dst,
                   TrilinosWrappers::MPI::BlockVector &tmp_distributed_disp)
  {
    if (parameters.nullspace_removal & NullspaceRemoval::angular_momentum)
      {
        // use_constant_density = false, remove net angular momentum
        remove_net_angular_momentum( false, relevant_dst, tmp_distributed_disp);
      }
    if (parameters.nullspace_removal & NullspaceRemoval::linear_momentum)
      {
        // use_constant_density = false, remove net momentum
        remove_net_linear_momentum( false, relevant_dst, tmp_distributed_disp);
      }
    if (parameters.nullspace_removal & NullspaceRemoval::net_rotation)
      {
        // use_constant_density = true, remove net rotation
        remove_net_angular_momentum( true, relevant_dst, tmp_distributed_disp);
      }
    if (parameters.nullspace_removal & NullspaceRemoval::net_translation)
      {
        // use_constant_density = true, remove net translation
        remove_net_linear_momentum( true, relevant_dst, tmp_distributed_disp);
      }
  }


  template <int dim>
  void 
  Simulator<dim>::
  remove_net_linear_momentum (const bool use_constant_density,
                              TrilinosWrappers::MPI::BlockVector &relevant_dst,
                              TrilinosWrappers::MPI::BlockVector &tmp_distributed_disp)
  {
    // compute and remove net linear momentum from displacement field, by computing
    // \int \rho (u + u_const) = 0

    QGauss<dim> quadrature(parameters.displacement_degree+1);
    const unsigned int n_q_points = quadrature.size();
    const unsigned int n_comp = introspection.n_compositional_fields;

    FEValues<dim> fe_values (*mapping, finite_element, quadrature,
                             update_values | update_quadrature_points | update_JxW_values);

    Tensor<1,dim> local_momentum;
    double local_mass = 0.0;

    // Vectors for evaluating the finite element solution
    std::vector<Tensor<1,dim> > local_displacements (n_q_points);

    const MaterialModel::MaterialProperties::Property requested_properties
      = MaterialModel::MaterialProperties::density;
    const MaterialModel::FieldDependences::Dependence field_dependences 
      = material_handler.get_field_dependences_for_evaluation(requested_properties);

    MaterialModel::MaterialModelInputs<dim> in(n_q_points, n_comp, 
                                               field_dependences, 
                                               requested_properties);
    MaterialModel::MaterialModelOutputs<dim> out(n_q_points, n_comp);

    // loop over all local cells
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
      {
        fe_values.reinit (cell);

        // get the displacement at each quadrature point
        fe_values[introspection.extractors.displacement].get_function_values (relevant_dst, 
                                                                              local_displacements);

        // get the density at each quadrature point if necessary
        if (!use_constant_density)
        {
          in.reinit (fe_values, qpd_handler, introspection, solution);
          material_handler.evaluate(in, out);
        }

        // actually compute the momentum and mass
        for (unsigned int q = 0; q < n_q_points; ++q)
        {
          // get the density at this quadrature point
          const double rho = (use_constant_density ? 1.0 : out.densities[q]);

          local_momentum += local_displacements[q] * rho * fe_values.JxW(q);
          local_mass += rho * fe_values.JxW(q);
        }
      }

    // Calculate the total mass and displacement correction
    const double mass = Utilities::MPI::sum( local_mass, mpi_communicator);
    Tensor<1,dim> displacement_correction = Utilities::MPI::sum(local_momentum, mpi_communicator)/mass;

    // We may only want to remove the nullspace for a single component, so zero out
    // the displacement correction if it is not selected by the NullspaceRemoval flag
    if (use_constant_density) // disable translation correction
      {
        if ( !(parameters.nullspace_removal & NullspaceRemoval::net_translation_x) )
          displacement_correction[0] = 0.0;  // don't correct x translation
        if ( !(parameters.nullspace_removal & NullspaceRemoval::net_translation_y) )
          displacement_correction[1] = 0.0;  // don't correct y translation
        if ( !(parameters.nullspace_removal & NullspaceRemoval::net_translation_z) && dim == 3 )
          displacement_correction[2] = 0.0;  // don't correct z translation
      }
    else // disable momentum correction
      {
        if ( !(parameters.nullspace_removal & NullspaceRemoval::linear_momentum_x) )
          displacement_correction[0] = 0.0;  // don't correct x translation
        if ( !(parameters.nullspace_removal & NullspaceRemoval::linear_momentum_y) )
          displacement_correction[1] = 0.0;  // don't correct y translation
        if ( !(parameters.nullspace_removal & NullspaceRemoval::linear_momentum_z) && dim == 3 )
          displacement_correction[2] = 0.0;  // don't correct z translation
      }

    // vector for storing the correction to the displacement field
    TrilinosWrappers::MPI::Vector correction(tmp_distributed_disp.block(introspection.block_indices.displacement));

    // Now construct a translation vector with the desired rate and subtract it from our vector
    internal::Translation<dim> translation( displacement_correction );
    interpolate_onto_displacement_system(translation, correction);
    tmp_distributed_disp.block(introspection.block_indices.displacement).add(-1.0,correction);

    // copy into the locally relevant vector
    relevant_dst.block(introspection.block_indices.displacement) =
      tmp_distributed_disp.block(introspection.block_indices.displacement);
  }


  template <int dim>
  RotationProperties<dim>
  Simulator<dim>::
  compute_net_angular_momentum(const bool use_constant_density,
                               const TrilinosWrappers::MPI::BlockVector &solution) const
  {
    // compute the momentum from displacement field, by computing
    // \int \rho u \cdot r_orth = \omega  * \int \rho x^2    ( 2 dimensions)
    // \int \rho r \times u =  I \cdot \omega  (3 dimensions)

    QGauss<dim> quadrature(parameters.displacement_degree+1);
    QGauss<dim-1> surface_quadrature(parameters.displacement_degree+1);

    const unsigned int n_q_points = quadrature.size();
    const unsigned int n_comp = introspection.n_compositional_fields;

    FEValues<dim> fe_values (*mapping, finite_element, quadrature,
                             update_values | update_quadrature_points | update_JxW_values);

    // moment of inertia and angular momentum for 3D
    SymmetricTensor<2,dim> local_moment_of_inertia;
    Tensor<1,dim> local_angular_momentum;
    std::vector<Tensor<1,dim> > local_displacements(n_q_points);

    // analogues to the moment of inertia and angular momentum for 2D
    double local_scalar_moment_of_inertia = 0.0;
    double local_scalar_angular_momentum = 0.0;

    // Structures for evaluating the displacements and the material model
    const MaterialModel::MaterialProperties::Property requested_properties
      = MaterialModel::MaterialProperties::density;
    const MaterialModel::FieldDependences::Dependence field_dependences 
      = material_handler.get_field_dependences_for_evaluation(requested_properties);

    MaterialModel::MaterialModelInputs<dim> in(n_q_points, n_comp, 
                                              field_dependences, 
                                              requested_properties);
    MaterialModel::MaterialModelOutputs<dim> out(n_q_points, n_comp);

    // loop over all local cells
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
      {
        fe_values.reinit (cell);

        if (!use_constant_density)
        {
          in.reinit (fe_values, qpd_handler, introspection, solution);
          material_handler.evaluate(in, out);
        }
        else
        {
          // Get the displacement at each quadrature point
          fe_values[introspection.extractors.displacement].get_function_values (solution, local_displacements);
        }

        // actually compute the moment of inertia and angular momentum
        for (unsigned int q = 0; q < n_q_points; ++q)
        {
          // get the position and density at this quadrature point
          const Point<dim> r_vec = fe_values.quadrature_point(q);
          const double rho = (use_constant_density ? 1.0 : out.densities[q]);
          const double JxW = fe_values.JxW(q);

          if (dim == 2)
          {
            // Get the displacement perpendicular to the position vector
            const Tensor<1,dim> r_perp = cross_product_2d(r_vec);

            local_scalar_angular_momentum += local_displacements[q] * r_perp * rho * JxW;
            local_scalar_moment_of_inertia += r_vec.norm_square() * rho * JxW;
          }
          else
          {
            const Tensor<1,dim> r_cross_v = cross_product_3d(r_vec, local_displacements[q]);

            local_angular_momentum += r_cross_v * rho * JxW;
            local_moment_of_inertia += (r_vec.norm_square() * unit_symmetric_tensor<dim>() - symmetrize(outer_product(r_vec,r_vec))) * rho * JxW;
          }
        }
      }

    RotationProperties<dim> properties;

    // Sum up the local contributions and solve for the overall rotation
    if (dim == 2)
    {
      properties.scalar_moment_of_inertia = Utilities::MPI::sum(local_scalar_moment_of_inertia, mpi_communicator);
      properties.scalar_angular_momentum = Utilities::MPI::sum(local_scalar_angular_momentum, mpi_communicator);
      properties.scalar_rotation = properties.scalar_angular_momentum / properties.scalar_moment_of_inertia;
    }
    else
    {
      properties.tensor_moment_of_inertia = Utilities::MPI::sum(local_moment_of_inertia,
                                                                mpi_communicator);
      properties.tensor_angular_momentum = Utilities::MPI::sum(local_angular_momentum, mpi_communicator );
      const SymmetricTensor<2,dim> inverse_moment (invert(Tensor<2,dim>(properties.tensor_moment_of_inertia)));
      properties.tensor_rotation = - inverse_moment * properties.tensor_angular_momentum;
    }

    return properties;
  }


  template <int dim>
  void 
  Simulator<dim>::
  remove_net_angular_momentum(const bool use_constant_density,
                              TrilinosWrappers::MPI::BlockVector &relevant_dst,
                              TrilinosWrappers::MPI::BlockVector &tmp_distributed_disp)
  {

    RotationProperties<dim> rotation_properties = compute_net_angular_momentum(use_constant_density,
                                                                               relevant_dst);

    // vector for storing the correction to the displacement field
    TrilinosWrappers::MPI::Vector correction(tmp_distributed_disp.block(introspection.block_indices.displacement));

    if (dim == 2)
      {
        // Now construct a rotation vector with the desired rate and subtract it from our vector
        const internal::Rotation<dim> rot(0);
        interpolate_onto_displacement_system(rot, correction);
        tmp_distributed_disp.block(introspection.block_indices.displacement).add(-1.0*rotation_properties.scalar_rotation,correction);
      }
    else
      {
        // Remove that rotation from the solution vector
        const internal::Rotation<dim> rot(rotation_properties.tensor_rotation);
        interpolate_onto_displacement_system(rot, correction);
        tmp_distributed_disp.block(introspection.block_indices.displacement).add(1.0,correction);
      }

    // copy into the locally relevant vector
    relevant_dst.block(introspection.block_indices.displacement) =
      tmp_distributed_disp.block(introspection.block_indices.displacement);
  }

}


// explicit instantiation of the functions we implement in this file
namespace elaspect
{
#define INSTANTIATE(dim) \
  template struct RotationProperties<dim>; \
  template void Simulator<dim>::remove_nullspace (TrilinosWrappers::MPI::BlockVector &, \
                                                  TrilinosWrappers::MPI::BlockVector &vector); \
  template void Simulator<dim>::setup_nullspace_constraints (AffineConstraints<double> &);

  ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
