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


#include <elaspect/simulator/assemblers/thermo_system.h>
#include <elaspect/boundary_heat_flux/interface.h>
#include <elaspect/mesh_deformation/handler.h>

namespace elaspect
{
  namespace Assemblers
  {
    namespace
    {
      /**
       * These functions implement a reduced form of the code from deal.II's TriaAccessor::measure().
       * In the 3d dG case, a call to face->measure() is not implemented for non-planar faces.
       * Since we only care about the scaling here, it is enough to have an approximation instead.
       * The 2d case remains unchanged.
       */
      double
      approximate_face_measure(const DoFHandler<2>::face_iterator &face)
      {
        return (face->vertex(0)-face->vertex(1)).norm();
      }

      double
      approximate_face_measure(const DoFHandler<3>::face_iterator &face)
      {
        const Tensor<1,3> v03 = face->vertex(3) - face->vertex(0);
        const Tensor<1,3> v12 = face->vertex(2) - face->vertex(1);
        const Tensor<1,3> twice_area = cross_product_3d(v03, v12);
        return 0.5 * twice_area.norm();
      }
    }


    template <int dim>
    void
    ThermoConvectionDiffusion<dim>::
    execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
            internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::ThermoSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::ThermoSystem<dim>&>(scratch_base);
      internal::Assembly::CopyData::ThermoSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::ThermoSystem<dim>&>(data_base);

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();

      const unsigned int T_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points = scratch.fe_values.n_quadrature_points;

      const bool use_ALE_method = this->get_parameters().use_ALE_method;

      const double time_step = this->get_timestep();

      for (unsigned int q = 0; q < n_q_points; ++q)
      {
        for (unsigned int i = 0, i_T = 0; i_T < T_dofs_per_cell; /*increment at end of loop*/)
        {
          if (fe.system_to_component_index(i).first == introspection.component_indices.temperature)
          {
            scratch.phi_T[i_T] = scratch.fe_values[introspection.extractors.temperature].value(i, q);
            scratch.grad_phi_T[i_T] = scratch.fe_values[introspection.extractors.temperature].gradient(i, q);
            ++i_T;
          }
          ++i;
        }

        const double density_c_P = scratch.material_model_outputs.densities[q] *
                                   scratch.material_model_outputs.specific_heat[q];

        const double latent_heat_LHS = scratch.heating_model_outputs.lhs_latent_heat_terms[q];

        AssertThrow(density_c_P + latent_heat_LHS >= 0,
                    ExcMessage ("The sum of density times c_P and the latent heat contribution "
                                "to the left hand side needs to be a non-negative quantity."));

        const double conductivity = scratch.material_model_outputs.thermal_conductivities[q];

        const double JxW = scratch.fe_values.JxW(q);

        for (unsigned int i = 0; i < T_dofs_per_cell; ++i)
        {
          data.local_rhs(i)
          += (scratch.phi_T[i]
              *
              (scratch.old_temperature_values[q] * (density_c_P + latent_heat_LHS)
               +
               time_step * scratch.heating_model_outputs.heating_source_terms[q])
             )
             * JxW;

          for (unsigned int j = 0; j < T_dofs_per_cell; ++j)
          {
            data.local_matrix(i,j)
            += ((scratch.phi_T[i] * scratch.phi_T[j]) * (density_c_P + latent_heat_LHS)
                +
                (time_step * conductivity * (scratch.grad_phi_T[i] * scratch.grad_phi_T[j]))
               ) 
               * JxW;
          }
        }

        if (use_ALE_method)
        {
          // add the convection term
          const Tensor<1,dim> beta = (scratch.displacement_increments[q] - 
                                      scratch.mesh_displacement_increments[q]) 
                                      * (density_c_P + latent_heat_LHS);

          for (unsigned int i = 0; i < T_dofs_per_cell; ++i)
            for (unsigned int j = 0; j < T_dofs_per_cell; ++j)
              data.local_matrix(i, j) += scratch.phi_T[i] * (beta * scratch.grad_phi_T[j]) * JxW;
        }
      }
    }


    template <int dim>
    MaterialModel::MaterialProperties::Property
    ThermoConvectionDiffusion<dim>::get_needed_material_properties() const
    {
      MaterialModel::MaterialProperties::Property
      properties = MaterialModel::MaterialProperties::density |
                   MaterialModel::MaterialProperties::specific_heat |
                   MaterialModel::MaterialProperties::thermal_conductivity;

      return properties;
    }


    template <int dim>
    void
    ThermoBoundaryFlux<dim>::
    execute(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base,
            internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      AssertThrow(this->get_parameters().use_ALE_method, ExcInternalError());

      internal::Assembly::Scratch::ThermoSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::ThermoSystem<dim>&>(scratch_base);
      internal::Assembly::CopyData::ThermoSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::ThermoSystem<dim>&>(data_base);
      
      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();

      const unsigned int face_no = scratch.face_number;
      const typename DoFHandler<dim>::face_iterator face = scratch.cell->face(face_no);

      const unsigned int n_face_q_points = scratch.fe_face_values->n_quadrature_points;
      const unsigned int T_dofs_per_cell = data.local_dof_indices.size();

      const unsigned int T_component = introspection.component_indices.temperature;
      const FEValuesExtractors::Scalar &T_extractor = introspection.extractors.temperature;

      const double time_step = this->get_timestep();

      if (this->get_fixed_heat_flux_boundary_indicators().find(face->boundary_id())
          != this->get_fixed_heat_flux_boundary_indicators().end())
      {
        for (unsigned int q = 0; q < n_face_q_points; ++q)
        {
          // Precompute the values of shape functions.
          for (unsigned int i = 0, i_T = 0; i_T < T_dofs_per_cell; /*increment at end of loop*/)
          {
            if (fe.system_to_component_index(i).first == T_component)
            {
              scratch.face_phi_T[i_T] = (*scratch.fe_face_values)[T_extractor].value(i, q);
              ++i_T;
            }
            ++i;
          }

          const Tensor<1,dim> heat_flux = 
            this->get_boundary_heat_flux().heat_flux(face->boundary_id(),
                                                     scratch.fe_face_values->quadrature_point(q),
                                                     scratch.fe_face_values->normal_vector(q));

          for (unsigned int i = 0; i < T_dofs_per_cell; ++i)
          {
            data.local_rhs(i)
            -= time_step * scratch.face_phi_T[i]
               * (heat_flux * scratch.fe_face_values->normal_vector(q))
               * scratch.fe_face_values->JxW(q);
          }
        }
      }
    }


    template <int dim>
    void
    ThermoSystemBoundaryFace<dim>::execute(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base,
                                           internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::ThermoSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::ThermoSystem<dim>&>(scratch_base);
      internal::Assembly::CopyData::ThermoSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::ThermoSystem<dim>&>(data_base);

      const Parameters<dim> &parameters = this->get_parameters();
      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();

      AssertThrow(parameters.use_ALE_method, ExcInternalError());

      const double time_step = this->get_timestep();

      const unsigned int face_no = scratch.face_number;
      const typename DoFHandler<dim>::face_iterator face = scratch.cell->face(face_no);

      const unsigned int n_q_points = scratch.fe_face_values->n_quadrature_points;
      const unsigned int T_dofs_per_cell = data.local_dof_indices.size();

      const unsigned int T_component = introspection.component_indices.temperature;
      const FEValuesExtractors::Scalar &T_extractor = introspection.extractors.temperature;

      if (this->get_fixed_temperature_boundary_indicators().find(face->boundary_id())
          != this->get_fixed_temperature_boundary_indicators().end())
      {
        // We are in the case of a Dirichlet temperature or composition boundary.
        // In the temperature case, impose the Dirichlet value weakly using a matrix term
        // and RHS term. In the composition case, Dirichlet conditions can only be imposed
        // on inflow boundaries, and we only have the flow-dependent terms, so we only
        // assemble the corresponding flow-dependent and matrix and RHS terms
        // if we are on an inflow boundary.
        for (unsigned int q = 0; q < n_q_points; ++q)
        {
          // precompute the values of shape functions and their gradients.
          for (unsigned int i = 0, i_T = 0; i_T < T_dofs_per_cell; /*increment at end of loop*/)
          {
            if (fe.system_to_component_index(i).first == T_component)
            {
              scratch.face_phi_T[i_T]      = (*scratch.fe_face_values)[T_extractor].value(i, q);
              scratch.face_grad_phi_T[i_T] = (*scratch.fe_face_values)[T_extractor].gradient(i, q);
              ++i_T;
            }
            ++i;
          }

          const double density_c_P = scratch.face_material_model_outputs.densities[q] *
                                     scratch.face_material_model_outputs.specific_heat[q];

          const double conductivity = scratch.face_material_model_outputs.thermal_conductivities[q];

          const double latent_heat_LHS = scratch.face_heating_model_outputs.lhs_latent_heat_terms[q];

          AssertThrow(density_c_P + latent_heat_LHS > 0,
                      ExcMessage("The sum of density times c_P and the latent heat contribution "
                                 "to the left hand side needs to be a positive quantity."));

          const double penalty = parameters.discontinuous_penalty
                                 * parameters.temperature_degree
                                 * parameters.temperature_degree
                                 / approximate_face_measure(face)
                                 * conductivity
                                 / (density_c_P + latent_heat_LHS);

          const double dirichlet_value = 
            this->get_boundary_temperature_manager().boundary_temperature(
              face->boundary_id(), scratch.fe_face_values->quadrature_point(q));

          // The discontinuous Galerkin method uses 2 types of jumps over edges:
          // undirected and directed jumps. Undirected jumps are dependent only
          // on the order of the numbering of cells. Directed jumps are dependent
          // on the direction of the flow. Thus the flow-dependent terms below are
          // only calculated if the edge is an inflow edge.
          const double du_n = (scratch.face_displacement_increments[q] -
                               scratch.face_mesh_displacement_increments[q])
                              * scratch.fe_face_values->normal_vector(q);
          const bool inflow = (du_n < 0.);

          for (unsigned int i = 0; i < T_dofs_per_cell; ++i)
          {
            data.local_rhs(i)
            += (- time_step * conductivity
                * scratch.face_grad_phi_T[i]
                * scratch.fe_face_values->normal_vector(q)
                * dirichlet_value

                + time_step
                * (density_c_P + latent_heat_LHS)
                * penalty
                * scratch.face_phi_T[i]
                * dirichlet_value

                + (inflow
                   ?
                   - (density_c_P + latent_heat_LHS)
                   * du_n
                   * dirichlet_value
                   * scratch.face_phi_T[i]
                   :
                   0.)
               )
               *
               scratch.fe_face_values->JxW(q);

            // local_matrix terms
            for (unsigned int j = 0; j < T_dofs_per_cell; ++j)
            {
              data.local_matrix(i,j)
              += (- time_step * conductivity
                  * scratch.face_grad_phi_T[i]
                  * scratch.fe_face_values->normal_vector(q)
                  * scratch.face_phi_T[j]

                  - time_step * conductivity
                  * scratch.face_grad_phi_T[j]
                  * scratch.fe_face_values->normal_vector(q)
                  * scratch.face_phi_T[i]

                  + time_step
                  * (density_c_P + latent_heat_LHS)
                  * penalty
                  * scratch.face_phi_T[i]
                  * scratch.face_phi_T[j]

                  + (inflow
                     ?
                     - (density_c_P + latent_heat_LHS)
                     * du_n
                     * scratch.face_phi_T[i]
                     * scratch.face_phi_T[j]
                     :
                     0.)
                 )
                 * scratch.fe_face_values->JxW(q);
            }
          }
        }
      }
      else
      {
        // Neumann temperature term - no non-zero contribution as only homogeneous
        // Neumann boundary conditions are implemented elsewhere for temperature
      }
    }


    template <int dim>
    MaterialModel::MaterialProperties::Property
    ThermoSystemBoundaryFace<dim>::get_needed_material_properties() const
    {
      MaterialModel::MaterialProperties::Property
      properties = MaterialModel::MaterialProperties::density |
                   MaterialModel::MaterialProperties::specific_heat |
                   MaterialModel::MaterialProperties::thermal_conductivity;

      return properties;
    }


    template <int dim>
    void
    ThermoSystemInteriorFace<dim>::execute(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base,
                                           internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::ThermoSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::ThermoSystem<dim>&>(scratch_base);
      internal::Assembly::CopyData::ThermoSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::ThermoSystem<dim>&>(data_base);

      const Parameters<dim> &parameters = this->get_parameters();
      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();

      AssertThrow(parameters.use_ALE_method, ExcInternalError());

      const double time_step = this->get_timestep();

      const typename DoFHandler<dim>::active_cell_iterator cell = scratch.cell;
      const unsigned int face_no = scratch.face_number;
      const typename DoFHandler<dim>::face_iterator face = cell->face(face_no);

      const unsigned int n_q_points = scratch.fe_face_values->n_quadrature_points;
      const unsigned int T_dofs_per_cell = data.local_dof_indices.size();

      const unsigned int T_component = introspection.component_indices.temperature;
      const FEValuesExtractors::Scalar &T_extractor = introspection.extractors.temperature;

      const typename DoFHandler<dim>::cell_iterator
      neighbor = cell->neighbor_or_periodic_neighbor(face_no);

      // note: "neighbor" defined above is NOT active_cell_iterator, so this includes cells that are refined
      // for example: cell with periodic boundary.
      Assert(neighbor.state() == IteratorState::valid, ExcInternalError());

      const bool cell_has_periodic_neighbor = cell->has_periodic_neighbor(face_no);

      if (!neighbor->has_children())
      {
        if (neighbor->level() == cell->level() &&
            neighbor->is_active() &&
            (((neighbor->is_locally_owned()) && (cell->index() < neighbor->index()))
             ||
             ((!neighbor->is_locally_owned()) && (cell->subdomain_id() < neighbor->subdomain_id()))))
        {
          // cell and neighbor are equal-sized, and cell has been chosen to assemble this face, 
          // so calculate from cell

          const unsigned int neighbor2 = 
            (cell_has_periodic_neighbor
             ?
             // how does the periodic neighbor talk about this cell ?
             cell->periodic_neighbor_of_periodic_neighbor(face_no)
             :
             // how does the neighbor talk about this cell?
             cell->neighbor_of_neighbor(face_no));

          // set up neighbor values
          scratch.neighbor_fe_face_values->reinit(neighbor, neighbor2);

          scratch.neighbor_face_material_model_inputs.reinit(*scratch.neighbor_fe_face_values,
                                                             this->get_qpd_handler(),
                                                             this->introspection(),
                                                             this->get_solution());

          this->get_material_handler().evaluate(scratch.neighbor_face_material_model_inputs,
                                                scratch.neighbor_face_material_model_outputs);

          this->get_heating_model_manager().evaluate(scratch.neighbor_face_material_model_inputs,
                                                     scratch.neighbor_face_material_model_outputs,
                                                     scratch.neighbor_face_heating_model_outputs);

          // get all dof indices on the neighbor, then extract those
          // that correspond to the temperature field
          neighbor->get_dof_indices(scratch.neighbor_dof_indices);
          for (unsigned int i = 0, i_T = 0; i_T < T_dofs_per_cell; /*increment at end of loop*/)
          {
            if (fe.system_to_component_index(i).first == T_component)
            {
              data.neighbor_dof_indices[face_no * GeometryInfo<dim>::max_children_per_face][i_T] = scratch.neighbor_dof_indices[i];
              ++i_T;
            }
            ++i;
          }
          data.assembled_matrices[face_no * GeometryInfo<dim>::max_children_per_face] = true;

          for (unsigned int q = 0; q < n_q_points; ++q)
          {
            // precompute the values of shape functions and their gradients.
            for (unsigned int i = 0, i_T = 0; i_T < T_dofs_per_cell; /*increment at end of loop*/)
            {
              if (fe.system_to_component_index(i).first == T_component)
              {
                scratch.face_phi_T[i_T]               = (*scratch.fe_face_values)[T_extractor].value(i, q);
                scratch.face_grad_phi_T[i_T]          = (*scratch.fe_face_values)[T_extractor].gradient(i, q);
                scratch.neighbor_face_phi_T[i_T]      = (*scratch.neighbor_fe_face_values)[T_extractor].value(i, q);
                scratch.neighbor_face_grad_phi_T[i_T] = (*scratch.neighbor_fe_face_values)[T_extractor].gradient(i, q);
                ++i_T;
              }
              ++i;
            }

            const double density_c_P = scratch.face_material_model_outputs.densities[q] *
                                       scratch.face_material_model_outputs.specific_heat[q];

            const double conductivity = scratch.face_material_model_outputs.thermal_conductivities[q];
            
            const double latent_heat_LHS = scratch.face_heating_model_outputs.lhs_latent_heat_terms[q];

            AssertThrow(density_c_P + latent_heat_LHS > 0,
                        ExcMessage("The sum of density times c_P and the latent heat contribution "
                                   "to the left hand side needs to be a positive quantity."));

            const double penalty = parameters.discontinuous_penalty
                                   * parameters.temperature_degree
                                   * parameters.temperature_degree
                                   / approximate_face_measure(face)
                                   * conductivity
                                   / (density_c_P + latent_heat_LHS);

            const double neighbor_density_c_P = 
              scratch.neighbor_face_material_model_outputs.densities[q] *
              scratch.neighbor_face_material_model_outputs.specific_heat[q];

            const double neighbor_conductivity = 
              scratch.neighbor_face_material_model_outputs.thermal_conductivities[q];

            const double neighbor_latent_heat_LHS =
              scratch.neighbor_face_heating_model_outputs.lhs_latent_heat_terms[q];

            AssertThrow(neighbor_density_c_P + neighbor_latent_heat_LHS > 0,
                        ExcMessage("The sum of density times c_P and the latent heat contribution "
                                   "to the left hand side on the neighbor needs to be a positive quantity."));

            const double neighbor_penalty = parameters.discontinuous_penalty
                                            * parameters.temperature_degree
                                            * parameters.temperature_degree
                                            / approximate_face_measure(neighbor->face(neighbor2))
                                            * neighbor_conductivity
                                            / (neighbor_density_c_P + neighbor_latent_heat_LHS);

            const double max_penalty = std::max(penalty, neighbor_penalty);

            const double max_density_c_P_and_latent_heat =
              std::max(density_c_P + latent_heat_LHS,
                       neighbor_density_c_P + neighbor_latent_heat_LHS);

            AssertThrow(numbers::is_finite(max_density_c_P_and_latent_heat),
                        ExcMessage("The maximum product of density and c_P plus latent LHS "
                                   "needs to be a finite quantity."));

             // The discontinuous Galerkin method uses 2 types of jumps over edges:
             // undirected and directed jumps. Undirected jumps are dependent only
             // on the order of the numbering of cells. Directed jumps are dependent
             // on the direction of the flow. Thus the flow-dependent terms below are
             // only calculated if the edge is an inflow edge.
            const double du_n = (scratch.face_displacement_increments[q] -
                                 scratch.face_mesh_displacement_increments[q])
                                * scratch.fe_face_values->normal_vector(q);
            const bool inflow = (du_n < 0.);

            for (unsigned int i = 0; i < T_dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < T_dofs_per_cell; ++j)
              {
                data.local_matrix(i, j)
                += (- 0.5 * time_step * conductivity
                    * scratch.face_grad_phi_T[i]
                    * scratch.fe_face_values->normal_vector(q)
                    * scratch.face_phi_T[j]

                    - 0.5 * time_step * conductivity
                    * scratch.face_grad_phi_T[j]
                    * scratch.fe_face_values->normal_vector(q)
                    * scratch.face_phi_T[i]

                    + time_step
                    * max_density_c_P_and_latent_heat
                    * max_penalty
                    * scratch.face_phi_T[i]
                    * scratch.face_phi_T[j]

                    - (inflow
                       ?
                       (density_c_P + latent_heat_LHS)
                       * du_n
                       * scratch.face_phi_T[i]
                       * scratch.face_phi_T[j]
                       :
                       0.)
                    )
                    * scratch.fe_face_values->JxW(q);
                
                data.local_matrices_int_ext[face_no * GeometryInfo<dim>::max_children_per_face](i, j)
                += (- 0.5 * time_step * neighbor_conductivity
                    * scratch.neighbor_face_grad_phi_T[j]
                    * scratch.fe_face_values->normal_vector(q)
                    * scratch.face_phi_T[i]

                    + 0.5 * time_step * conductivity
                    * scratch.face_grad_phi_T[i]
                    * scratch.fe_face_values->normal_vector(q)
                    * scratch.neighbor_face_phi_T[j]

                    - time_step
                    * max_density_c_P_and_latent_heat
                    * max_penalty
                    * scratch.neighbor_face_phi_T[j]
                    * scratch.face_phi_T[i]

                    + (inflow
                       ?
                       (density_c_P + latent_heat_LHS)
                       * du_n
                       * scratch.face_phi_T[i]
                       * scratch.neighbor_face_phi_T[j]
                       :
                       0.)
                   )
                   * scratch.fe_face_values->JxW(q);

                data.local_matrices_ext_int[face_no * GeometryInfo<dim>::max_children_per_face](i, j)
                += (+ 0.5 * time_step * conductivity
                    * scratch.face_grad_phi_T[j]
                    * scratch.fe_face_values->normal_vector(q)
                    * scratch.neighbor_face_phi_T[i]

                    - 0.5 * time_step * neighbor_conductivity
                    * scratch.neighbor_face_grad_phi_T[i]
                    * scratch.fe_face_values->normal_vector(q)
                    * scratch.face_phi_T[j]

                    - time_step
                    * max_density_c_P_and_latent_heat
                    * max_penalty
                    * scratch.face_phi_T[j]
                    * scratch.neighbor_face_phi_T[i]

                    - (!inflow
                       ?
                       (neighbor_density_c_P + neighbor_latent_heat_LHS)
                       * du_n
                       * scratch.neighbor_face_phi_T[i]
                       * scratch.face_phi_T[j]
                       :
                       0.)
                   )
                   * scratch.fe_face_values->JxW(q);

                data.local_matrices_ext_ext[face_no * GeometryInfo<dim>::max_children_per_face](i, j)
                += (+ 0.5 * time_step * neighbor_conductivity
                    * scratch.neighbor_face_grad_phi_T[i]
                    * scratch.fe_face_values->normal_vector(q)
                    * scratch.neighbor_face_phi_T[j]

                    + 0.5 * time_step * neighbor_conductivity
                    * scratch.neighbor_face_grad_phi_T[j]
                    * scratch.fe_face_values->normal_vector(q)
                    * scratch.neighbor_face_phi_T[i]

                    + time_step
                    * max_density_c_P_and_latent_heat
                    * max_penalty
                    * scratch.neighbor_face_phi_T[i]
                    * scratch.neighbor_face_phi_T[j]

                    + (!inflow
                       ?
                       (neighbor_density_c_P + neighbor_latent_heat_LHS)
                       * du_n
                       * scratch.neighbor_face_phi_T[i]
                       * scratch.neighbor_face_phi_T[j]
                       :
                       0.)
                   )
                   * scratch.fe_face_values->JxW(q);
              }
            }
          }
        }
        else
        {
          // neighbor is taking responsibility for assembly of this face, because
          // either (1) neighbor is coarser, or
          //        (2) neighbor is equally-sized and
          //           (a) neighbor is on a different subdomain, with lower subdmain_id(), or
          //           (b) neighbor is on the same subdomain and has lower index().
        }
      }
      // neighbor has children, so always assemble from here.
      else
      {
        const unsigned int neighbor2 =
          (cell_has_periodic_neighbor
           ?
           cell->periodic_neighbor_face_no(face_no)
           :
           cell->neighbor_face_no(face_no));

        // Loop over subfaces. We know that the neighbor is finer, so we could loop over the subfaces of the current
        // face. but if we are at a periodic boundary, then the face of the current cell has no children, so instead use
        // the children of the periodic neighbor's corresponding face since we know that the letter does indeed have
        // children (because we know that the neighbor is refined).
        typename DoFHandler<dim>::face_iterator neighbor_face = neighbor->face(neighbor2);
        for (unsigned int subface_no = 0; subface_no < neighbor_face->n_children(); ++subface_no)
        {
          const typename DoFHandler<dim>::active_cell_iterator neighbor_child
            = (cell_has_periodic_neighbor
               ?
               cell->periodic_neighbor_child_on_subface(face_no, subface_no)
               :
               cell->neighbor_child_on_subface(face_no, subface_no));

          // set up subface values
          scratch.fe_subface_values->reinit(cell, face_no, subface_no);

          (*scratch.fe_subface_values)[introspection.extractors.displacement].get_function_values(
            this->get_solution(), scratch.face_displacement_increments);
          (*scratch.fe_subface_values)[introspection.extractors.displacement].get_function_values(
            this->get_mesh_deformation_handler().get_mesh_displacement_increments(), scratch.face_mesh_displacement_increments);

          scratch.face_material_model_inputs.reinit(*scratch.fe_subface_values,
                                                    this->get_qpd_handler(),
                                                    this->introspection(),
                                                    this->get_solution());

          this->get_material_handler().evaluate(scratch.face_material_model_inputs,
                                                scratch.face_material_model_outputs);

          this->get_heating_model_manager().evaluate(scratch.face_material_model_inputs,
                                                     scratch.face_material_model_outputs,
                                                     scratch.face_heating_model_outputs);

          // set up neighbor values
          scratch.neighbor_fe_face_values->reinit(neighbor_child, neighbor2);

          scratch.neighbor_face_material_model_inputs.reinit(*scratch.neighbor_fe_face_values,
                                                             this->get_qpd_handler(),
                                                             this->introspection(),
                                                             this->get_solution());

          this->get_material_handler().evaluate(scratch.neighbor_face_material_model_inputs,
                                                scratch.neighbor_face_material_model_outputs);

          this->get_heating_model_manager().evaluate(scratch.neighbor_face_material_model_inputs,
                                                     scratch.neighbor_face_material_model_outputs,
                                                     scratch.neighbor_face_heating_model_outputs);

          // get all dof indices on the neighbor, then extract those
          // that correspond to temperature field.
          neighbor_child->get_dof_indices(scratch.neighbor_dof_indices);
          for (unsigned int i = 0, i_T = 0; i_T < T_dofs_per_cell; /*increment at end of loop*/)
          {
            if (fe.system_to_component_index(i).first == T_component)
            {
              data.neighbor_dof_indices[face_no * GeometryInfo<dim>::max_children_per_face + subface_no][i_T] = scratch.neighbor_dof_indices[i];
              ++i_T;
            }
            ++i;
          }
          data.assembled_matrices[face_no * GeometryInfo<dim>::max_children_per_face + subface_no] = true;

          for (unsigned int q = 0; q < n_q_points; ++q)
          {
            // precompute the values of shape functions and their gradients
            for (unsigned int i = 0, i_T = 0; i_T < T_dofs_per_cell; /*increment at end of loop*/)
            {
              if (fe.system_to_component_index(i).first == T_component)
              {
                scratch.face_phi_T[i_T]               = (*scratch.fe_subface_values)[T_extractor].value(i, q);
                scratch.face_grad_phi_T[i_T]          = (*scratch.fe_subface_values)[T_extractor].gradient(i, q);
                scratch.neighbor_face_phi_T[i_T]      = (*scratch.neighbor_fe_face_values)[T_extractor].value(i, q);
                scratch.neighbor_face_grad_phi_T[i_T] = (*scratch.neighbor_fe_face_values)[T_extractor].gradient(i, q);
                ++i_T;
              }
              ++i;
            }

            const double density_c_P = scratch.face_material_model_outputs.densities[q] *
                                       scratch.face_material_model_outputs.specific_heat[q];

            const double conductivity = scratch.face_material_model_outputs.thermal_conductivities[q];

            const double latent_heat_LHS = scratch.face_heating_model_outputs.lhs_latent_heat_terms[q];

            AssertThrow(density_c_P + latent_heat_LHS > 0,
                        ExcMessage("The sum of density times c_P and the latent heat contribution "
                                   "to the left hand side needs to be a positive quantity."));

            const double penalty = parameters.discontinuous_penalty
                                   * parameters.temperature_degree
                                   * parameters.temperature_degree
                                   / approximate_face_measure(face)
                                   * conductivity
                                   / (density_c_P + latent_heat_LHS);

            const double neighbor_density_c_P =
              scratch.neighbor_face_material_model_outputs.densities[q] *
              scratch.neighbor_face_material_model_outputs.specific_heat[q];

            const double neighbor_conductivity =
              scratch.neighbor_face_material_model_outputs.thermal_conductivities[q];

            const double neighbor_latent_heat_LHS =
              scratch.neighbor_face_heating_model_outputs.lhs_latent_heat_terms[q];

            AssertThrow(neighbor_density_c_P + neighbor_latent_heat_LHS > 0,
                        ExcMessage("The sum of density times c_P and the latent heat contribution "
                                   "to the left hand side on the neighbor needs to be a positive quantity."));

            const double neighbor_penalty = parameters.discontinuous_penalty
                                            * parameters.temperature_degree
                                            * parameters.temperature_degree
                                            / approximate_face_measure(neighbor_child->face(neighbor2))
                                            * neighbor_conductivity
                                            / (neighbor_density_c_P + neighbor_latent_heat_LHS);

            const double max_penalty = std::max(penalty, neighbor_penalty);

            const double max_density_c_P_and_latent_heat =
              std::max(density_c_P + latent_heat_LHS,
                       neighbor_density_c_P + neighbor_latent_heat_LHS);

            AssertThrow(numbers::is_finite(max_density_c_P_and_latent_heat),
                        ExcMessage("The maximum product of density and c_P plus latent heat LHS "
                                   "on the neighbor needs to be a finite quantity."));
            
            // The discontinuous Galerkin method uses 2 types of jumps over edges:
            // undirected and directed jumps. Undirected jumps are dependent only
            // on the order of the numbering of cells. Directed jumps are dependent
            // on the direction of the flow. Thus the flow-dependent terms below are
            // only calculated if the edge is an inflow edge.
            const double du_n = (scratch.face_displacement_increments[q] -
                                 scratch.face_mesh_displacement_increments[q])
                                * scratch.fe_subface_values->normal_vector(q);
            const bool inflow = (du_n < 0.);

            for (unsigned int i = 0; i < T_dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < T_dofs_per_cell; ++j)
              {
                data.local_matrix(i, j)
                += (- 0.5 * time_step * conductivity
                    * scratch.face_grad_phi_T[i]
                    * scratch.fe_subface_values->normal_vector(q)
                    * scratch.face_phi_T[j]

                    - 0.5 * time_step * conductivity
                    * scratch.face_grad_phi_T[j]
                    * scratch.fe_subface_values->normal_vector(q)
                    * scratch.face_phi_T[i]

                    + time_step
                    * max_density_c_P_and_latent_heat
                    * max_penalty
                    * scratch.face_phi_T[i]
                    * scratch.face_phi_T[j]

                    - (inflow
                       ?
                       (density_c_P + latent_heat_LHS)
                       * du_n
                       * scratch.face_phi_T[i]
                       * scratch.face_phi_T[j]
                       :
                       0.)
                   )
                   * scratch.fe_subface_values->JxW(q);

                data.local_matrices_int_ext[face_no * GeometryInfo<dim>::max_children_per_face + subface_no](i, j)
                += (- 0.5 * time_step * neighbor_conductivity
                    * scratch.neighbor_face_grad_phi_T[j]
                    * scratch.fe_subface_values->normal_vector(q)
                    * scratch.face_phi_T[i]

                    + 0.5 * time_step * conductivity
                    * scratch.face_grad_phi_T[i]
                    * scratch.fe_subface_values->normal_vector(q)
                    * scratch.neighbor_face_phi_T[j]

                    - time_step
                    * max_density_c_P_and_latent_heat
                    * max_penalty
                    * scratch.neighbor_face_phi_T[j]
                    * scratch.face_phi_T[i]

                    + (inflow
                       ?
                       (density_c_P + latent_heat_LHS)
                       * du_n
                       * scratch.face_phi_T[i]
                       * scratch.neighbor_face_phi_T[j]
:
                       0.)
            )
                   * scratch.fe_subface_values->JxW(q);

                data.local_matrices_ext_int[face_no * GeometryInfo<dim>::max_children_per_face + subface_no](i, j)
                += (+ 0.5 * time_step * conductivity
                    * scratch.face_grad_phi_T[j]
                    * scratch.fe_subface_values->normal_vector(q)
* scratch.neighbor_face_phi_T[i]

                    - 0.5 * time_step * neighbor_conductivity
                    * scratch.neighbor_face_grad_phi_T[i]
                    * scratch.fe_subface_values->normal_vector(q)
                    * scratch.face_phi_T[j]

                    - time_step
                    * max_density_c_P_and_latent_heat
                    * max_penalty
* scratch.face_phi_T[j]
                    * scratch.neighbor_face_phi_T[i]

                    - (!inflow
                       ?
                       (neighbor_density_c_P + neighbor_latent_heat_LHS)
                       * du_n
                       * scratch.neighbor_face_phi_T[i]
                       * scratch.face_phi_T[j]
                       :
                       0.)
                   )
                   * scratch.fe_subface_values->JxW(q);

                data.local_matrices_ext_ext[face_no * GeometryInfo<dim>::max_children_per_face + subface_no](i,j)
                += (+ 0.5 * time_step * neighbor_conductivity
                    * scratch.neighbor_face_grad_phi_T[i]
                    * scratch.fe_subface_values->normal_vector(q)
                    * scratch.neighbor_face_phi_T[j]

                    + 0.5 * time_step * neighbor_conductivity
                    * scratch.neighbor_face_grad_phi_T[j]
                    * scratch.fe_subface_values->normal_vector(q)
                    * scratch.neighbor_face_phi_T[i]

                    + time_step
                    * max_density_c_P_and_latent_heat
                    * max_penalty
                    * scratch.neighbor_face_phi_T[i]
                    * scratch.neighbor_face_phi_T[j]

                    + (!inflow
                       ?
                       (neighbor_density_c_P + neighbor_latent_heat_LHS)
                       * du_n
                       * scratch.neighbor_face_phi_T[i]
                       * scratch.neighbor_face_phi_T[j]
                       :
                       0.)
                   )
                   * scratch.fe_subface_values->JxW(q);
              }
            }
          }
        }
      }
    }


    template <int dim>
    MaterialModel::MaterialProperties::Property
    ThermoSystemInteriorFace<dim>::get_needed_material_properties() const
    {
      MaterialModel::MaterialProperties::Property
      properties = MaterialModel::MaterialProperties::density |
                   MaterialModel::MaterialProperties::specific_heat |
                   MaterialModel::MaterialProperties::thermal_conductivity;

      return properties;
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace Assemblers
  {
#define INSTANTIATE(dim) \
    template class ThermoConvectionDiffusion<dim>; \
    template class ThermoBoundaryFlux<dim>; \
    template class ThermoSystemBoundaryFace<dim>; \
    template class ThermoSystemInteriorFace<dim>;

    ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
