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


#include <elaspect/time_stepping/conduction_time_step.h>
#include <elaspect/material_model/handler.h>

#include <deal.II/base/quadrature_lib.h>

namespace elaspect
{
  namespace TimeStepping
  {
    template <int dim>
    double
    ConductionTimeStep<dim>::get_next_time_step_size () const
    {
      double min_local_conduction_timestep = std::numeric_limits<double>::max();

      const QIterated<dim> quadrature_formula (QTrapezoid<1>(),
                                               this->get_parameters().displacement_degree);

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values |
                               update_gradients |
                               update_quadrature_points);

      const unsigned int n_q_points = quadrature_formula.size();

      const MaterialModel::MaterialProperties::Property requested_properties 
        = MaterialModel::MaterialProperties::density |
          MaterialModel::MaterialProperties::specific_heat |
          MaterialModel::MaterialProperties::thermal_conductivity;

      const MaterialModel::FieldDependences::Dependence field_dependences 
        = this->get_material_handler().get_field_dependences_for_evaluation (requested_properties);

      MaterialModel::MaterialModelInputs<dim> in (n_q_points,
                                                  this->introspection().n_compositional_fields,
                                                  field_dependences,
                                                  requested_properties);
      MaterialModel::MaterialModelOutputs<dim> out (n_q_points,
                                                    this->introspection().n_compositional_fields);

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          in.reinit (fe_values,
                     this->get_qpd_handler(),
                     this->introspection(),
                     this->get_solution());

          this->get_material_handler().evaluate (in, out);

          for (unsigned int q = 0; q < n_q_points; ++q)
          {
            const double k   = out.thermal_conductivities[q];
            const double rho = out.densities[q];
            const double c_p = out.specific_heat[q];

            Assert (rho * c_p > 0,
                    ExcMessage ("The product of density and c_P needs to be a "
                                "non-negative quantity."));

            const double thermal_diffusivity = k / (rho * c_p);

            if (thermal_diffusivity > 0)
            {
              min_local_conduction_timestep = std::min (min_local_conduction_timestep,
                                                        this->get_parameters().CFL_number * pow(cell->minimum_vertex_distance(), 2.)
                                                        / thermal_diffusivity);
            }
          }
        }

      const double min_conduction_timestep = Utilities::MPI::min (min_local_conduction_timestep, this->get_mpi_communicator());

      AssertThrow (min_conduction_timestep > 0,
                   ExcMessage("The time step length for the each time step needs to be positive, "
                              "but the computed step length was: " + std::to_string(min_conduction_timestep) + ". "
                              "Please check for non-positive material properties."));

      return min_conduction_timestep;
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace TimeStepping
  {
    ELASPECT_REGISTER_TIME_STEPPING_MODEL(ConductionTimeStep,
                                      "conduction time step",
                                      "This model computes the conduction time step as the minimum "
                                      "over all cells of $ CFL h^2 \\cdot \\rho C_p / k$, "
                                      "where k is the thermal conductivity. This plugin will always "
                                      "request advancing to the next time step.")
  }
}
