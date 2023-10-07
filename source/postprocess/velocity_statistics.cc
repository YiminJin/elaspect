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


#include <elaspect/postprocess/velocity_statistics.h>

namespace elaspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    VelocityStatistics<dim>::execute(TableHandler &statistics)
    {
      FEValues<dim> fe_values(this->get_mapping(),
                              this->get_fe(),
                              this->get_quadrature_formula(),
                              update_values | update_JxW_values);

      const unsigned int n_q_points = fe_values.n_quadrature_points;
      std::vector<Tensor<1,dim>> displacement_increments(n_q_points);

      double local_du2_integral = 0;
      double local_max_du = 0;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values[this->introspection().extractors.displacement].get_function_values(
            this->get_solution(), displacement_increments);

          for (unsigned int q = 0; q < n_q_points; ++q)
          {
            const double du2 = displacement_increments[q] * displacement_increments[q];
            local_du2_integral += du2 * fe_values.JxW(q);
            local_max_du = std::max(du2, local_max_du);
          }
        }

      const double global_du2_integral
        = Utilities::MPI::sum(local_du2_integral, this->get_mpi_communicator());
      const double global_max_du
        = Utilities::MPI::max(local_max_du, this->get_mpi_communicator());

      const double rms = std::sqrt(global_du2_integral / this->get_volume());
      const double one_over_dt = 1. / ( this->get_timestep() *
                                        ( this->convert_output_to_years() ? 
                                          year_in_seconds : 1.0 ) );

      const std::string units = (this->convert_output_to_years() ? "m/year" : "m/s");
      const std::vector<std::string> column_names = { "RMS velocity (" + units + ")",
                                                      "Max. velocity (" + units + ")"
                                                    };

      statistics.add_value(column_names[0], rms * one_over_dt);
      statistics.add_value(column_names[1], global_max_du * one_over_dt);

      // also make sure that the other columns filled by this object
      // all show up with sufficient accuracy and in scientific notation
      for (auto &column : column_names)
      {
        statistics.set_precision (column, 8);
        statistics.set_scientific (column, true);
      }

      std::ostringstream output;
      output.precision(3);
      output << rms * one_over_dt << ", "
             << global_max_du * one_over_dt;

      return std::pair<std::string, std::string>("RMS, max velocity (" + units + "):",
                                                 output.str());
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace Postprocess
  {
    ELASPECT_REGISTER_POSTPROCESSOR(VelocityStatistics,
                                "velocity statistics",
                                "A postprocessor that computes some statistics about the "
                                "velocity field.")
  }
}
