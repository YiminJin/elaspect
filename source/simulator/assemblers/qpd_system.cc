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


#include <elaspect/simulator/assemblers/qpd_system.h>
#include <elaspect/simulator.h>

namespace elaspect
{
  namespace Assemblers
  {
    template <int dim>
    void
    QPDProjection<dim>::execute(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base,
                                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      AssertThrow(!this->get_parameters().use_ALE_method, ExcInternalError());

      internal::Assembly::Scratch::QPDSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::QPDSystem<dim>&>(scratch_base);
      internal::Assembly::CopyData::QPDSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::QPDSystem<dim>&>(data_base);

      const FiniteElement<dim> &fe = this->get_fe();

      const unsigned int n_q_points = scratch.fe_values.n_quadrature_points;
      const unsigned int field_dofs_per_cell = data.local_dof_indices[0].size();

      // get the values of shape functions through the first QPD component,
      // since all compositional/physical fields share one base element
      const unsigned int field0_component = scratch.qpd_fields[0].component_index;
      const FEValuesExtractors::Scalar field0_extractor = scratch.qpd_fields[0].scalar_extractor;

      for (unsigned int q = 0; q < n_q_points; ++q)
      {
        // precompute the values of shape functions
        for (unsigned int i = 0, i_field = 0; i_field < field_dofs_per_cell; /*increment at end of loop*/)
        {
          if (fe.system_to_component_index(i).first == field0_component)
          {
            scratch.phi_field[i_field] = scratch.fe_values[field0_extractor].value(i, q);
            ++i_field;
          }
          ++i;
        }

        // do the actual assembly.
        for (unsigned int i = 0; i < field_dofs_per_cell; ++i)
        {
          for (unsigned int j = 0; j < field_dofs_per_cell; ++j)
            data.local_matrix(i, j)
            += (scratch.phi_field[i] * scratch.phi_field[j])
               * scratch.fe_values.JxW(q);

          for (unsigned int field = 0; field < scratch.qpd_fields.size(); ++field)
            data.local_rhs[field](i)
            += (scratch.phi_field[i] * scratch.old_field_values[field][q])
               * scratch.fe_values.JxW(q);
        }
      }
    }


    template <int dim>
    void
    QPDAdvection<dim>::execute(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base,
                               internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      AssertThrow(this->get_parameters().use_ALE_method, ExcInternalError());

      internal::Assembly::Scratch::QPDSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::QPDSystem<dim>&>(scratch_base);
      internal::Assembly::CopyData::QPDSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::QPDSystem<dim>&>(data_base);

      const FiniteElement<dim> &fe = this->get_fe();

      const unsigned int n_q_points = scratch.fe_values.n_quadrature_points;
      const unsigned int field_dofs_per_cell = data.local_dof_indices[0].size();

      // get the values and gradients of shape functions through the first QPD component,
      // since all compositional/physical fields share one base element
      const unsigned int field0_component = scratch.qpd_fields[0].component_index;
      const FEValuesExtractors::Scalar &field0_extractor = scratch.qpd_fields[0].scalar_extractor;

      for (unsigned int q = 0; q < n_q_points; ++q)
      {
        // precompute the values and gradients of shape functions
        for (unsigned int i = 0, i_field = 0; i_field < field_dofs_per_cell; /*increment at end of loop*/)
        {
          if (fe.system_to_component_index(i).first == field0_component)
          {
            scratch.grad_phi_field[i_field] = scratch.fe_values[field0_extractor].gradient(i, q);
            scratch.phi_field[i_field]      = scratch.fe_values[field0_extractor].value(i, q);
            ++i_field;
          }
          ++i;
        }

        // compute the relative displacement increment
        const Tensor<1,dim> du = scratch.displacement_increments[q] -
                                 scratch.mesh_displacement_increments[q];

        // do the actual assembly.
        for (unsigned int i = 0; i < field_dofs_per_cell; ++i)
        {
          for (unsigned int j = 0; j < field_dofs_per_cell; ++j)
            data.local_matrix(i, j)
            += (scratch.phi_field[i] * 
                (scratch.phi_field[j] + du * scratch.grad_phi_field[j]))
               * scratch.fe_values.JxW(q);

          for (unsigned int field = 0; field < scratch.qpd_fields.size(); ++field)
            data.local_rhs[field](i)
            += (scratch.phi_field[i] * scratch.old_field_values[field][q])
               * scratch.fe_values.JxW(q);
        }
      }
    }


    template <int dim>
    void
    QPDAdvectionBoundaryFace<dim>::execute(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base,
                                           internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      AssertThrow(this->get_parameters().use_ALE_method, ExcInternalError());

      internal::Assembly::Scratch::QPDSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::QPDSystem<dim>&>(scratch_base);
      internal::Assembly::CopyData::QPDSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::QPDSystem<dim>&>(data_base);

      // Dirichlet boundary conditions can only be applied to compositional fields
      if (!scratch.qpd_fields[0].is_compositional_field())
        return;

      const FiniteElement<dim> &fe = this->get_fe();

      const unsigned int face_no = scratch.face_number;
      const typename DoFHandler<dim>::face_iterator face = scratch.cell->face(face_no);

      const unsigned int n_q_points = scratch.fe_face_values->n_quadrature_points;
      const unsigned int field_dofs_per_cell = data.local_dof_indices[0].size();

      // get the values of shape functions through the first QPD component,
      // since all compositional/physical fields share one base element
      const unsigned int field0_component = scratch.qpd_fields[0].component_index;
      const FEValuesExtractors::Scalar &field0_extractor = scratch.qpd_fields[0].scalar_extractor;

      if (this->get_fixed_composition_boundary_indicators().find(face->boundary_id()) 
          != this->get_fixed_composition_boundary_indicators().end())
      {
        // We are in the case of a Dirichlet composition boundary.
        for (unsigned int q = 0; q < n_q_points; ++q)
        {
          // Precompute the values of shape functions and their gradients.
          for (unsigned int i = 0, i_field = 0; i_field < field_dofs_per_cell; /*increment at end of loop*/)
          {
            if (fe.system_to_component_index(i).first == field0_component)
            {
              scratch.face_phi_field[i_field] = (*scratch.fe_face_values)[field0_extractor].value(i, q);
              ++i_field;
            }
            ++i;
          }
          
          // Dirichlet conditions can only be imposed on inflow boundaries, so we only 
          // assemble the corresponding matrix and RHS terms if we are on an inflow 
          // boundary.
          const double du_n = (scratch.face_displacement_increments[q] -
                               scratch.face_mesh_displacement_increments[q])
                              * scratch.fe_face_values->normal_vector(q);

          if (du_n < 0)
          {
            for (unsigned int field = 0; field < scratch.qpd_fields.size(); ++field)
            {
              const double dirichlet_value = 
                this->get_boundary_composition_manager().boundary_composition(
                    face->boundary_id(),
                    scratch.fe_face_values->quadrature_point(q),
                    scratch.qpd_fields[field].compositional_variable);

              for (unsigned int i = 0; i < field_dofs_per_cell; ++i)
                data.local_rhs[field](i) 
                -= dirichlet_value * scratch.face_phi_field[i]
                   * du_n * scratch.fe_face_values->JxW(q);
            }

            for (unsigned int i = 0; i < field_dofs_per_cell; ++i)
              for (unsigned int j = 0; j < field_dofs_per_cell; ++j)
                data.local_matrix(i, j)
                -= scratch.face_phi_field[i] * scratch.face_phi_field[j]
                   * du_n * scratch.fe_face_values->JxW(q);
          }
        }
      }
      else
      {
        // Newmann boundary -- no non-zero contribution as only homogeneous Neumann boundary
        // conditions are implemented
      }
    }


    template <int dim>
    void
    QPDAdvectionInteriorFace<dim>::execute(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base,
                                           internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      AssertThrow(this->get_parameters().use_ALE_method, ExcInternalError());

      internal::Assembly::Scratch::QPDSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::QPDSystem<dim>&>(scratch_base);
      internal::Assembly::CopyData::QPDSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::QPDSystem<dim>&>(data_base);

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();

      const typename DoFHandler<dim>::active_cell_iterator cell = scratch.cell;
      const unsigned int face_no = scratch.face_number;

      const unsigned int n_q_points = scratch.fe_face_values->n_quadrature_points;
      const unsigned int field_dofs_per_cell = data.local_dof_indices[0].size();

      // get the values of shape functions through the first QPD component,
      // since all compositional/physical fields share one base element
      const unsigned int field0_component = scratch.qpd_fields[0].component_index;
      const FEValuesExtractors::Scalar &field0_extractor = scratch.qpd_fields[0].scalar_extractor;

      const typename DoFHandler<dim>::cell_iterator
      neighbor = cell->neighbor_or_periodic_neighbor(face_no);
      Assert(neighbor.state() == IteratorState::valid,
             ExcInternalError());
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
            (cell->has_periodic_neighbor(face_no)
             ?
             // how does the periodic neighbor talk about this cell?
             cell->periodic_neighbor_of_periodic_neighbor(face_no)
             :
             // how does the neighbor talk about this cell?
             cell->neighbor_of_neighbor(face_no));

          // set up neighbor values
          scratch.neighbor_fe_face_values->reinit(neighbor, neighbor2);

          // get all dof indices on the neighbor, then extract those 
          // that correspond to the solution_field we are interested in
          neighbor->get_dof_indices(scratch.neighbor_dof_indices);
          for (unsigned int i = 0, i_field = 0; i_field < field_dofs_per_cell; /*increment at end of loop*/)
          {
            if (fe.system_to_component_index(i).first == field0_component)
            {
              data.neighbor_dof_indices[face_no * GeometryInfo<dim>::max_children_per_face][i_field] = scratch.neighbor_dof_indices[i];
              ++i_field;
            }
            ++i;
          }
          data.assembled_matrices[face_no * GeometryInfo<dim>::max_children_per_face] = true;

          for (unsigned int q = 0; q < n_q_points; ++q)
          {
            // precompute the values of shape functions.
            for (unsigned int i = 0, i_field = 0; i_field < field_dofs_per_cell; /*increment at end of loop*/)
            {
              if (fe.system_to_component_index(i).first == field0_component)
              {
                scratch.face_phi_field[i_field]          = (*scratch.fe_face_values)[field0_extractor].value(i, q);
                scratch.neighbor_face_phi_field[i_field] = (*scratch.neighbor_fe_face_values)[field0_extractor].value(i, q);
                ++i_field;
              }
              ++i;
            }

            const double du_n = (scratch.face_displacement_increments[q] - 
                                 scratch.face_mesh_displacement_increments[q])
                                * scratch.fe_face_values->normal_vector(q);

            // The discontinuous Galerkin method uses 2 types of jumps over edges:
            // undirected and directed jumps. Undirected jumps are dependent only
            // on the order of the numbering of cells. Directed jumps are dependent
            // on the direction of the flow. Thus the flow-dependent terms below are
            // only calculated if the edge is an inflow edge.
            for (unsigned int i = 0; i < field_dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < field_dofs_per_cell; ++j)
              {
                if (du_n < 0)
                {
                  data.local_matrix(i, j)
                  -= (scratch.face_phi_field[i] * scratch.face_phi_field[j])
                     * du_n * scratch.fe_face_values->JxW(q);
                
                  data.local_matrices_int_ext[face_no * GeometryInfo<dim>::max_children_per_face](i, j)
                  += (scratch.face_phi_field[i] * scratch.neighbor_face_phi_field[j])
                     * du_n * scratch.fe_face_values->JxW(q);
                }
                else
                {
                  data.local_matrices_ext_int[face_no * GeometryInfo<dim>::max_children_per_face](i, j)
                  -= (scratch.neighbor_face_phi_field[i] * scratch.face_phi_field[j])
                     * du_n * scratch.fe_face_values->JxW(q);

                  data.local_matrices_ext_ext[face_no * GeometryInfo<dim>::max_children_per_face](i, j)
                  += (scratch.neighbor_face_phi_field[i] * scratch.neighbor_face_phi_field[j])
                     * du_n * scratch.fe_face_values->JxW(q);
                }
              }
            }
          }
        }
        else
        {
          // neighbor is taking responsibility for assembly of this face, because
          // either (1) neighbor is coarser, or
          //        (2) neighbor is equally-sized and
          //            (a) neighbor is on a different subdomain, with lower subdomain_id(), or
          //            (b) neighbor is on the same subdomain and has lower index().
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
        // face. But if we are at a periodic boundary, then the face of the current cell has no children, so instead use
        // the children of the periodic neighbor's corresponding face since we know that the later does indeed have
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
            this->get_solution(),
            scratch.face_displacement_increments);
          (*scratch.fe_subface_values)[introspection.extractors.displacement].get_function_values(
            this->get_mesh_deformation_handler().get_mesh_displacement_increments(), 
            scratch.face_mesh_displacement_increments);

          // set up neighbor values
          scratch.neighbor_fe_face_values->reinit(neighbor_child, neighbor2);

          // get all dof indices on the neighbor, then extract those
          // that correspond to the solution_field we are interested in
          neighbor_child->get_dof_indices(scratch.neighbor_dof_indices);
          for (unsigned int i = 0, i_field = 0; i_field < field_dofs_per_cell; /*increment at end of loop*/)
          {
            if (fe.system_to_component_index(i).first == field0_component)
            {
              data.neighbor_dof_indices[face_no * GeometryInfo<dim>::max_children_per_face + subface_no][i_field] = scratch.neighbor_dof_indices[i];
              ++i_field;
            }
            ++i;
          }

          data.assembled_matrices[face_no * GeometryInfo<dim>::max_children_per_face + subface_no] = true;

          for (unsigned int q = 0; q < n_q_points; ++q)
          {
            // precompute the values of shape functions.
            for (unsigned int i = 0, i_field = 0; i_field < field_dofs_per_cell; /*increment at end of loop*/)
            {
              if (fe.system_to_component_index(i).first == field0_component)
              {
                scratch.face_phi_field[i_field]          = (*scratch.fe_subface_values)[field0_extractor].value(i, q);
                scratch.neighbor_face_phi_field[i_field] = (*scratch.neighbor_fe_face_values)[field0_extractor].value(i, q);
                ++i_field;
              }
              ++i;
            }
            
            const double du_n = ( scratch.face_displacement_increments[q] -
                                  scratch.face_mesh_displacement_increments[q] )
                                * scratch.fe_subface_values->normal_vector(q);

            // The discontinuous Galerkin method uses 2 types of jumps over edges:
            // undirected and directed jumps. Undirected jumps are dependent only
            // on the order of the numbering of cells. Directed jumps are dependent
            // on the direction of the flow. Thus the flow-dependent terms below are
            // only calculated if the edge is an inflow edge.
            for (unsigned int i = 0; i < field_dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < field_dofs_per_cell; ++j)
              {
                if (du_n < 0)
                {
                  data.local_matrix(i, j)
                  -= (scratch.face_phi_field[i] * scratch.face_phi_field[j])
                     * du_n * scratch.fe_subface_values->JxW(q);

                  data.local_matrices_int_ext[face_no * GeometryInfo<dim>::max_children_per_face + subface_no](i, j)
                  += (scratch.face_phi_field[i] * scratch.neighbor_face_phi_field[j])
                     * du_n * scratch.fe_subface_values->JxW(q);
                }
                else
                {
                  data.local_matrices_ext_int[face_no * GeometryInfo<dim>::max_children_per_face + subface_no](i, j)
                  -= (scratch.neighbor_face_phi_field[i] * scratch.face_phi_field[j])
                     * du_n * scratch.fe_subface_values->JxW(q);

                  data.local_matrices_ext_ext[face_no * GeometryInfo<dim>::max_children_per_face + subface_no](i, j)
                  += (scratch.neighbor_face_phi_field[i] * scratch.neighbor_face_phi_field[j])
                     * du_n * scratch.fe_subface_values->JxW(q);
                }
              }
            }
          }
        }
      }
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace Assemblers
  {
#define INSTANTIATE(dim) \
    template class QPDProjection<dim>; \
    template class QPDAdvection<dim>; \
    template class QPDAdvectionBoundaryFace<dim>; \
    template class QPDAdvectionInteriorFace<dim>;

    ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
