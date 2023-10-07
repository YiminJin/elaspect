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


#ifndef _elaspect_material_model_average_h
#define _elaspect_material_model_average_h

#include <elaspect/global.h>
#include <elaspect/material_model/io_interface.h>

#include <deal.II/base/quadrature.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/lac/full_matrix.h>

namespace elaspect
{
  namespace MaterialModel
  {
    namespace MaterialAveraging
    {
      using namespace dealii;

      enum AveragingOperation
      {
        none,
        arithmetic_average,
        harmonic_average,
        geometric_average,
        pick_largest,
        project_to_Q1,
        log_average
      };

      std::string get_averaging_operation_names ();

      AveragingOperation parse_averaging_operation_name (const std::string &s);

      template <int dim>
      void average (const AveragingOperation operation,
                    const typename DoFHandler<dim>::active_cell_iterator &cell,
                    const Quadrature<dim> &quadrature_formula,
                    const Mapping<dim> &mapping,
                    const bool average_tangent_moduli,
                    MaterialModelOutputs<dim> &values_out);

      void average_property (const AveragingOperation operation,
                             const FullMatrix<double> &projection_matrix,
                             const FullMatrix<double> &expansion_matrix,
                             std::vector<double>      &values_out);

      template <int dim>
      void compute_projection_matrix (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                      const Quadrature<dim> &quadrature_formula,
                                      const Mapping<dim>    &mapping,
                                      FullMatrix<double>    &projection_matrix,
                                      FullMatrix<double>    &expansion_matrix);
    }
  }
}

#endif
