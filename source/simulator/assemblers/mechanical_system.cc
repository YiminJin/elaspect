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


#include <elaspect/simulator/assemblers/mechanical_system.h>
#include <elaspect/gravity_model/interface.h>
#include <elaspect/boundary_traction/interface.h>
#include <elaspect/material_model/utilities.h>

#include <deal.II/physics/elasticity/standard_tensors.h>

namespace elaspect
{
  namespace Assemblers
  {
    template <int dim>
    void
    MechanicalQuasiStaticTerms<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::MechanicalSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::MechanicalSystem<dim>&>(scratch_base);
      internal::Assembly::CopyData::MechanicalSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::MechanicalSystem<dim>&>(data_base);

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int u_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points  = scratch.fe_values.n_quadrature_points;

      for (unsigned int q = 0; q < n_q_points; ++q)
      {
        for (unsigned int i = 0, i_u = 0; i_u < u_dofs_per_cell; /*increment at end of loop*/)
        {
          if (fe.system_to_component_index(i).first < dim)
          {
            scratch.phi_u[i_u] = scratch.fe_values[introspection.extractors.displacement].value(i, q);
            scratch.symgrad_phi_u[i_u] = scratch.fe_values[introspection.extractors.displacement].symmetric_gradient(i, q);
            ++i_u;           
          }
          ++i;
        }

        const Tensor<1,dim> gravity = 
          this->get_gravity_model().gravity_vector (scratch.fe_values.quadrature_point(q));
        const double density = scratch.material_model_outputs.densities[q];
        const double JxW = scratch.fe_values.JxW(q);

        SymmetricTensor<2,dim> stress_for_rhs;
        if ((this->get_nonlinear_iteration() > 0) &&
            (this->constitutive_relation() & ConstitutiveRelation::plasticity))
        {
          stress_for_rhs = scratch.material_model_inputs.qpd_cell->get_symmetric_tensor(q, introspection.qpd_indicators.current_stress);
        }
        else
        {
          stress_for_rhs = scratch.material_model_inputs.old_stress[q];

          if (this->get_timestep_number() > 0)
          {
            if (this->constitutive_relation() & ConstitutiveRelation::viscosity)
            {
              const double eta = scratch.material_model_outputs.viscosities[q];
              const double Gve = MaterialUtilities::get_shear_modulus(scratch.material_model_outputs.tangent_moduli[q]);

              const double p_old = trace(stress_for_rhs) / dim;
              SymmetricTensor<2,dim> s_old = deviator(stress_for_rhs) * ((eta - Gve * this->get_timestep()) / eta);

              stress_for_rhs = s_old + p_old * unit_symmetric_tensor<dim>();
            }

            if (this->constitutive_relation() & ConstitutiveRelation::thermal_expansion)
            {
              const double T_old = scratch.material_model_inputs.qpd_cell->get_scalar(q, introspection.qpd_indicators.old_temperature);
              const double alpha_dT = scratch.material_model_outputs.thermal_expansivities[q]
                                * (scratch.material_model_inputs.temperature[q] - T_old);

              stress_for_rhs -= scratch.material_model_outputs.tangent_moduli[q]
                                * (alpha_dT * unit_symmetric_tensor<dim>());
            }
          }
        }

        for (unsigned int i = 0; i < u_dofs_per_cell; ++i)
        {
          data.local_rhs(i) += ( density * (gravity * scratch.phi_u[i]) -
                                 stress_for_rhs * scratch.symgrad_phi_u[i] 
                               ) * JxW;
          
          if (this->reassemble_tangent_matrix())
          {
            for (unsigned int j = 0; j < u_dofs_per_cell; ++j)
              data.local_matrix(i, j) += ( scratch.symgrad_phi_u[i] *
                                           scratch.material_model_outputs.tangent_moduli[q] *
                                           scratch.symgrad_phi_u[j]
                                         ) * JxW;
          }
        }
      }
    }


    template <int dim>
    MaterialModel::FieldDependences::Dependence
    MechanicalQuasiStaticTerms<dim>::get_field_dependences() const
    {
      MaterialModel::FieldDependences::Dependence
      dependences = MaterialModel::FieldDependences::none;

      if (!(this->constitutive_relation() & ConstitutiveRelation::plasticity))
      {
        dependences |= MaterialModel::FieldDependences::old_stress;

        if (this->constitutive_relation() & ConstitutiveRelation::thermal_expansion)
          dependences |= MaterialModel::FieldDependences::temperature;
      }

      return dependences;
    }


    template <int dim>
    MaterialModel::MaterialProperties::Property
    MechanicalQuasiStaticTerms<dim>::get_needed_material_properties() const
    {
      MaterialModel::MaterialProperties::Property
      properties = MaterialModel::MaterialProperties::density;

      if (this->reassemble_tangent_matrix())
        properties |= MaterialModel::MaterialProperties::tangent_modulus;

      if (this->constitutive_relation() & ConstitutiveRelation::viscosity)
        properties |= MaterialModel::MaterialProperties::viscosity;

      if (this->constitutive_relation() & ConstitutiveRelation::thermal_expansion)
        properties |= MaterialModel::MaterialProperties::thermal_expansivity;

      return properties;
    }


    template <int dim>
    void
    MechanicalBoundaryTraction<dim>::execute (internal::Assembly::Scratch::ScratchBase<dim> &scratch_base,
                                         internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::MechanicalSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::MechanicalSystem<dim>&>(scratch_base);
      internal::Assembly::CopyData::MechanicalSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::MechanicalSystem<dim>&>(data_base);

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();

      const unsigned int u_dofs_per_cell = data.local_dof_indices.size();

      // see if any of the faces are traction boundaries for which
      // we need to assemble force terms for the right hand side
      const typename DoFHandler<dim>::face_iterator face = scratch.cell->face(scratch.face_number);

      if (this->get_boundary_traction().find(face->boundary_id())
          != this->get_boundary_traction().end())
      {
        for (unsigned int q = 0; q < scratch.fe_face_values.n_quadrature_points; ++q)
        {
          const Tensor<1,dim> traction = this->get_boundary_traction()
            .find(face->boundary_id())->second->boundary_traction (face->boundary_id(),
                                                                   scratch.fe_face_values.quadrature_point(q),
                                                                   scratch.fe_face_values.normal_vector(q));

          for (unsigned int i = 0, i_u = 0; i_u < u_dofs_per_cell; /*increment at end of loop*/)
          {
            if (fe.system_to_component_index(i).first < dim)
            {
              data.local_rhs(i_u) += scratch.fe_face_values[introspection.extractors.displacement].value(i,q) *
                                     traction *
                                     scratch.fe_face_values.JxW(q);
              ++i_u;
            }
            ++i;
          }
        }
      }
    }
  } // namespace Assemblers
} // namespace elaspect


// explicit instantiations
namespace elaspect
{
  namespace Assemblers
  {
#define INSTANTIATE(dim) \
    template class MechanicalQuasiStaticTerms<dim>; \
    template class MechanicalBoundaryTraction<dim>;

    ELASPECT_INSTANTIATE(INSTANTIATE)
  }
}
