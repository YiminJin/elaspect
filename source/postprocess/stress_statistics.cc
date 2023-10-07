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


#include <elaspect/postprocess/stress_statistics.h>

namespace elaspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    StressStatistics<dim>::execute(TableHandler &statistics)
    {
      FEValues<dim> fe_values(this->get_mapping(),
                              this->get_fe(),
                              this->get_quadrature_formula(),
                              update_JxW_values);

      const unsigned int n_q_points = fe_values.n_quadrature_points;
      const unsigned int stress_indicator = this->introspection().qpd_indicators.old_stress;

      // compute the integral quantities by quadrature
      std::vector<double> local_stress_integrals(SymmetricTensor<2,dim>::n_independent_components);

      typename DoFHandler<dim>::active_cell_iterator
      dof_cell = this->get_dof_handler().begin_active(),
      dof_endc = this->get_dof_handler().end();
      typename QPDHandler<dim>::active_cell_iterator
      qpd_cell = this->get_qpd_handler().begin_active();
      for (; dof_cell != dof_endc; ++dof_cell, ++qpd_cell)
        if (dof_cell->is_locally_owned())
        {
          fe_values.reinit(dof_cell);

          for (unsigned int q = 0; q < n_q_points; ++q)
          {
            SymmetricTensor<2,dim> stress = qpd_cell->get_symmetric_tensor(q, stress_indicator);
            for (unsigned int c = 0; c < SymmetricTensor<2,dim>::n_independent_components; ++c)
              local_stress_integrals[c] += stress.access_raw_entry(c) * fe_values.JxW(q);
          }
        }

      // compute the global average for each stress component
      std::vector<double> global_stress_integrals(SymmetricTensor<2,dim>::n_independent_components),
                          global_avg_stress(SymmetricTensor<2,dim>::n_independent_components);
      Utilities::MPI::sum(local_stress_integrals,
                          this->get_mpi_communicator(),
                          global_stress_integrals);
      for (unsigned int c = 0; c < SymmetricTensor<2,dim>::n_independent_components; ++c)
        global_avg_stress[c] = global_stress_integrals[c] / this->get_volume();

      // compute min/max by simply looping over the elements of the
      // solution vector
      std::vector<double> local_min_stress(SymmetricTensor<2,dim>::n_independent_components,
                                           std::numeric_limits<double>::max());
      std::vector<double> local_max_stress(SymmetricTensor<2,dim>::n_independent_components,
                                           std::numeric_limits<double>::lowest());

      for (unsigned int c = 0; c < SymmetricTensor<2,dim>::n_independent_components; ++c)
      {
        const unsigned int block_idx = this->introspection().block_indices.stress[c];
        IndexSet range = this->get_solution().block(block_idx).locally_owned_elements();
        for (unsigned int i = 0; i < range.n_elements(); ++i)
        {
          const unsigned int idx = range.nth_index_in_set(i);
          const double val = this->get_solution().block(block_idx)(idx);

          local_min_stress[c] = std::min<double>(local_min_stress[c], val);
          local_max_stress[c] = std::max<double>(local_max_stress[c], val);
        }
      }

      // do the reductions over all processors
      std::vector<double> global_min_stress(SymmetricTensor<2,dim>::n_independent_components,
                                            std::numeric_limits<double>::max());
      std::vector<double> global_max_stress(SymmetricTensor<2,dim>::n_independent_components,
                                            std::numeric_limits<double>::lowest());

      Utilities::MPI::min(local_min_stress,
                          this->get_mpi_communicator(),
                          global_min_stress);
      Utilities::MPI::max(local_max_stress,
                          this->get_mpi_communicator(),
                          global_max_stress);

      // finally produce something for the statistics file
      std::vector<std::string> stress_component_names;
      stress_component_names.emplace_back("stress_xx (Pa)");
      stress_component_names.emplace_back("stress_yy (Pa)");
      if (dim == 2)
        stress_component_names.emplace_back("stress_xy (Pa)");
      else
      {
        stress_component_names.emplace_back("stress_zz (Pa)");
        stress_component_names.emplace_back("stress_xy (Pa)");
        stress_component_names.emplace_back("stress_xz (Pa)");
        stress_component_names.emplace_back("stress_yz (Pa)");
      }
      for (unsigned int c = 0; c < SymmetricTensor<2,dim>::n_independent_components; ++c)
      {
        statistics.add_value("Minimal " + stress_component_names[c],
                             global_min_stress[c]);
        statistics.add_value("Average " + stress_component_names[c],
                             global_avg_stress[c]);
        statistics.add_value("Maximal " + stress_component_names[c],
                             global_max_stress[c]);

        // make sure that the columns filled by this object all show up
        // with sufficient accuracy and in scientific notation
        const std::string columns[] = { "Minimal " + stress_component_names[c],
                                        "Average " + stress_component_names[c],
                                        "Maximal " + stress_component_names[c] };
        for (const auto &col : columns)
        {
          statistics.set_precision(col, 8);
          statistics.set_scientific(col, true);
        }
      }

      std::ostringstream output;
      output.precision(4);
      for (unsigned int c = 0; c < SymmetricTensor<2,dim>::n_independent_components; ++c)
      {
        output << global_min_stress[c] << '/'
               << global_avg_stress[c] << '/'
               << global_max_stress[c];
        if (c+1 != SymmetricTensor<2,dim>::n_independent_components)
          output << " // ";
      }

      return std::pair<std::string, std::string>("Stress min/avg/max (Pa):",
                                                 output.str());
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace Postprocess
  {
    ELASPECT_REGISTER_POSTPROCESSOR(StressStatistics,
                                "stress statistics",
                                "A postprocessor that computes some statistics about "
                                "the stress field. In particular, it computes maximal and "
                                "minimal values of each stress component, as well as the "
                                "total mass as defined by the integral "
                                "$m_i(t) = \\int_\\Omega c_i(\\mathbf x,t) \\; \\text{d}x$.")
  }
}
