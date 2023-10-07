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


#include <elaspect/postprocess/temperature_statistics.h>

namespace elaspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    TemperatureStatistics<dim>::execute(TableHandler &statistics)
{
      FEValues<dim> fe_values(this->get_mapping(),
                              this->get_fe(),
                              this->get_quadrature_formula(),
                              update_values | update_JxW_values);

      const unsigned int n_q_points = fe_values.n_quadrature_points;
      std::vector<double> temperature_values(n_q_points);

      double local_temperature_integral = 0;

      // compute the integral quantities by quadrature
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values[this->introspection().extractors.temperature].get_function_values(
            this->get_solution(), temperature_values);

          for (unsigned int q = 0; q < n_q_points; ++q)
            local_temperature_integral += temperature_values[q] * fe_values.JxW(q);
        }

      // compute min/max by simply looping over the elements of the
      // solution vector. the reason is that minimum and maximum are
      // usually attained at the  boundary, and so taking their
      // values at Gauss quadrature points gives an inaccurate
      // picture of their true values
      double local_min_temperature = std::numeric_limits<double>::max();
      double local_max_temperature = std::numeric_limits<double>::lowest();
      const unsigned int temperature_block = this->introspection().block_indices.temperature;
      IndexSet range = this->get_solution().block(temperature_block).locally_owned_elements();
      for (unsigned int i = 0; i < range.n_elements(); ++i)
      {
        const unsigned int idx = range.nth_index_in_set(i);
        const double val = this->get_solution().block(temperature_block)(idx);

        local_min_temperature = std::min<double>(local_min_temperature, val);
        local_max_temperature = std::max<double>(local_max_temperature, val);
      }

      const double global_temperature_integral =
        Utilities::MPI::sum(local_temperature_integral, this->get_mpi_communicator());
      const double global_mean_temperature = global_temperature_integral / this->get_volume();

      double global_min_temperature = 0;
      double global_max_temperature = 0;
      // do the reductions that are min/max operations. do them in
      // one communication by multiplying one value by -1
      {
        double local_values[2] = { -local_min_temperature, local_max_temperature };
        double global_values[2];

        Utilities::MPI::max(local_values, this->get_mpi_communicator(), global_values);

        global_min_temperature = -global_values[0];
        global_max_temperature = global_values[1];
      }

      statistics.add_value("Minimal temperature (K)", global_min_temperature);
      statistics.add_value("Average temperature (K)", global_mean_temperature);
      statistics.add_value("Maximal temperature (K)", global_max_temperature);

      // also make sure that the other columns filled by this object
      // all show up with sufficient accuracy and in scientific notation
      {
        const char *columns[] = { "Minimal temperature (K)",
                                  "Average temperature (K)",
                                  "Maximal temperature (K)"
                                };
        for (auto &column : columns)
        {
          statistics.set_precision(column, 8);
          statistics.set_scientific(column, true);
        }    
      }

      std::ostringstream output;
      output.precision(4);
      output << global_min_temperature << '/'
             << global_temperature_integral / this->get_volume() << '/'
             << global_max_temperature;

      return std::pair<std::string, std::string> ("Temperature min/avg/max (K):",
                                                  output.str());   
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace Postprocess
  {
    ELASPECT_REGISTER_POSTPROCESSOR(TemperatureStatistics,
                                "temperature statistics",
                                "A postprocessor that computes some statistics about "
                                "the temperature field.")
  }
}
