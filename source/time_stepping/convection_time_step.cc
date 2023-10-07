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


#include <elaspect/time_stepping/convection_time_step.h>

#include <deal.II/base/quadrature_lib.h>

namespace elaspect
{
  namespace TimeStepping
  {
    template <int dim>
    double
    ConvectionTimeStep<dim>::get_next_time_step_size () const
    {
      if (this->get_timestep_number() == 0)
        return std::numeric_limits<double>::max();

      const QIterated<dim> quadrature_formula (QTrapezoid<1>(),
                                               this->get_parameters().displacement_degree);

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values);

      const unsigned int n_q_points = quadrature_formula.size();

      std::vector<Tensor<1,dim> > displacement_values (n_q_points);
      std::vector<Tensor<1,dim> > old_displacement_values (n_q_points);

      double max_local_speed_over_meshsize = 0;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          fe_values[this->introspection().extractors.displacement].get_function_values (this->get_solution(),
                                                                                         displacement_values);
          fe_values[this->introspection().extractors.displacement].get_function_values (this->get_old_solution(),
                                                                                         old_displacement_values);

          double max_local_velocity = 0;
          for (unsigned int q = 0; q < n_q_points; ++q)
          {
            const Tensor<1,dim> velocity = (displacement_values[q] - old_displacement_values[q])
                                           / this->get_timestep();
            max_local_velocity = std::max (max_local_velocity, velocity.norm());
          }

          max_local_speed_over_meshsize = std::max (max_local_speed_over_meshsize,
                                                    max_local_velocity / cell->minimum_vertex_distance());
        }

      const double max_global_speed_over_meshsize =
        Utilities::MPI::max (max_local_speed_over_meshsize, this->get_mpi_communicator());

      double min_convection_timestep = std::numeric_limits<double>::max();

      if (max_global_speed_over_meshsize != 0.0)
        min_convection_timestep = this->get_parameters().CFL_number / (this->get_parameters().temperature_degree * max_global_speed_over_meshsize);

      AssertThrow (min_convection_timestep > 0,
                   ExcMessage("The time step length for the each time step needs to be positive, "
                              "but the computed step length was: " + std::to_string(min_convection_timestep) + ". "
                              "Please check for non-positive material properties."));

      return min_convection_timestep;
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace TimeStepping
  {
    ELASPECT_REGISTER_TIME_STEPPING_MODEL(ConvectionTimeStep,
                                      "convection time step",
                                      "This model computes the convection time step as "
                                      "$ CFL / \\max \\| u \\| / h$ over all cells, "
                                      "where $u$ is the velocity and $h$ is the product of mesh size "
                                      "and temperature polynomial degree.")
  }
}
