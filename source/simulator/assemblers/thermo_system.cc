#include <elaspect/simulator/assemblers/thermo_system.h>
#include <elaspect/boundary_heat_flux/interface.h>

namespace elaspect
{
  namespace Assemblers
  {
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
      const unsigned int fe_order = fe.base_element(introspection.base_elements.temperature).degree;

      const bool use_ALE_method = this->get_parameters().use_ALE_method;

      const double time_step = this->get_timestep();

      // Compute the SUPG parameter tau defined in "On Discontinuity-Capturing 
      // Methods for Convection-Diffusion Equations" by Volker John and Petr
      // Knobloch. 
      // delta_k = h / (2 \|u\| p) * (coth(Pe) - 1/Pe)
      // Pe = \| u \| h/(2 p eps)
      double tau = 0.0;
      if (use_ALE_method)
      {
        double norm_of_advection_term = 0.0;
        double max_conductivity_on_cell = 0.0;

        for (unsigned int q = 0; q < n_q_points; ++q)
        {
          const Tensor<1,dim> relative_displacement_increment = 
            scratch.displacement_increments[q] - scratch.mesh_displacement_increments[q];

          norm_of_advection_term = 
            std::max(relative_displacement_increment.norm() / time_step *
                     (scratch.material_model_outputs.densities[q] * 
                      scratch.material_model_outputs.specific_heat[q] +
                      scratch.heating_model_outputs.lhs_latent_heat_terms[q]),
                     norm_of_advection_term);

          max_conductivity_on_cell = 
            std::max(scratch.material_model_outputs.thermal_conductivities[q],
                     max_conductivity_on_cell);
        }

        const double h = scratch.cell->diameter();
        const double eps = max_conductivity_on_cell;

        const double peclet_times_eps = norm_of_advection_term * h / (2.0 * fe_order);

        // Instead of Pe < 1, we check Pe * eps < eps as eps can be 0:
        if (!(peclet_times_eps == 0.0 || peclet_times_eps < eps))
        {
          // To avoid a division by zero, increase eps slightly. The actual value is not
          // important, as long as the result is still a valid number. Note that this
          // is only important if \|u\| and eps are zero.
          const double peclet = peclet_times_eps / (eps + 1e-100);
          const double coth_peclet = (1.0 + exp(-2.0 * peclet)) / (1.0 - exp(-2.0 * peclet));

          tau = h / (2.0 * norm_of_advection_term * fe_order) * (coth_peclet - 1.0 / peclet);

          Assert(tau >= 0, ExcMessage("tau for SUPG needs to be a nonnegative constant."));
        }
      }

      for (unsigned int q = 0; q < n_q_points; ++q)
      {
        for (unsigned int i = 0, i_T = 0; i_T < T_dofs_per_cell; /*increment at end of loop*/)
        {
          if (fe.system_to_component_index(i).first == introspection.component_indices.temperature)
          {
            scratch.phi_T[i_T] = scratch.fe_values[introspection.extractors.temperature].value(i, q);
            scratch.grad_phi_T[i_T] = scratch.fe_values[introspection.extractors.temperature].gradient(i, q);
            if (use_ALE_method)
              scratch.laplacian_phi_T[i_T] = trace(scratch.fe_values[introspection.extractors.temperature].hessian(i, q));
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
          += ((scratch.phi_T[i] * scratch.old_temperature_values[q]
               * (density_c_P + latent_heat_LHS))
               + 
               (time_step * scratch.phi_T[i] 
                * scratch.heating_model_outputs.heating_source_terms[q])
             )
             * JxW;
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
          const Tensor<1,dim> beta = (scratch.displacement_increments[q] - 
                                      scratch.mesh_displacement_increments[q]) 
                                      * (density_c_P + latent_heat_LHS) 
                                      / time_step;

          for (unsigned int i = 0; i < T_dofs_per_cell; ++i)
          {
            data.local_rhs(i)
            += (tau * (beta * scratch.grad_phi_T[i])
                * 
                (scratch.old_temperature_values[q] * (density_c_P + latent_heat_LHS) 
                 + 
                 time_step * scratch.heating_model_outputs.heating_source_terms[q])
               )
               * JxW;

            for (unsigned int j = 0; j < T_dofs_per_cell; ++j)
              data.local_matrix(i, j)
              += (time_step * scratch.phi_T[i] * (beta * scratch.grad_phi_T[j])
                  +
                  tau * (beta * scratch.grad_phi_T[i]) 
                  *
                  (scratch.phi_T[j] * (density_c_P + latent_heat_LHS)
                   +
                   time_step * (beta * scratch.grad_phi_T[j] -
                                conductivity * scratch.laplacian_phi_T[j]))
                 )
                 * JxW;
          }
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
      internal::Assembly::Scratch::ThermoSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::ThermoSystem<dim>&>(scratch_base);
      internal::Assembly::CopyData::ThermoSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::ThermoSystem<dim>&>(data_base);
      
      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();

      const unsigned int face_no = scratch.face_number;
      const typename DoFHandler<dim>::face_iterator face = scratch.cell->face(face_no);

      const unsigned int n_face_q_points = scratch.fe_face_values.n_quadrature_points;
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
              scratch.face_phi_T[i_T] = scratch.fe_face_values[T_extractor].value(i, q);
              ++i_T;
            }
            ++i;
          }

          const Tensor<1,dim> heat_flux = 
            this->get_boundary_heat_flux().heat_flux(face->boundary_id(),
                                                     scratch.fe_face_values.quadrature_point(q),
                                                     scratch.fe_face_values.normal_vector(q));

          for (unsigned int i = 0; i < T_dofs_per_cell; ++i)
          {
            data.local_rhs(i)
            -= time_step * scratch.face_phi_T[i]
               * (heat_flux * scratch.fe_face_values.normal_vector(q))
               * scratch.fe_face_values.JxW(q);
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
    template class ThermoConvectionDiffusion<dim>; \
    template class ThermoBoundaryFlux<dim>;

    ELASPECT_INSTANTIATE(INSTANTIATE)
  }
}
