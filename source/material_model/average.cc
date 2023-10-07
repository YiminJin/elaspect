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


#include <elaspect/material_model/average.h>

#include <deal.II/fe/fe_q.h>

#ifdef DEBUG
#ifdef ASPECT_USE_FP_EXCEPTIONS
#include <cfenv>
#endif
#endif

namespace elaspect
{
  namespace MaterialModel
  {
    namespace MaterialAveraging
    {
      std::string get_averaging_operation_names ()
      {
        return "none|arithmetic average|harmonic average|geometric average|pick largest|project to Q1|log average";
      }


      AveragingOperation parse_averaging_operation_name (const std::string &s)
      {
        if (s == "none")
          return none;
        else if (s == "arithmetic average")
          return arithmetic_average;
        else if (s == "harmonic average")
          return harmonic_average;
        else if (s == "geometric average")
          return geometric_average;
        else if (s == "pick largest")
          return pick_largest;
        else if (s == "project to Q1")
          return project_to_Q1;
        else if (s == "log average")
          return log_average;
        else
          AssertThrow (false,
                       ExcMessage ("The value <" + s + "> for a material "
                                   "averaging operation is not one of the "
                                   "valid values."));

        return none;
      }


      void average_property (const AveragingOperation operation,
                             const FullMatrix<double> &projection_matrix,
                             const FullMatrix<double> &expansion_matrix,
                             std::vector<double>      &values_out)
      {
        // if an output field has not been filled (because it was
        // not requested), then simply do nothing -- no harm no foul
        if (values_out.size() == 0)
          return;

#ifdef DEBUG
#ifdef ASPECT_USE_FP_EXCEPTIONS
        // disable floating point exceptions while averaging. Errors will be reported
        // as soon as somebody will try to use the averaged values later.
        fedisableexcept(FE_DIVBYZERO|FE_INVALID);
#endif
#endif

        const unsigned int N = values_out.size();
        const unsigned int P = expansion_matrix.n();
        Assert ((P==0) || (/*dim=2*/ P==4) || (/*dim=3*/ P==8),
                ExcInternalError());
        Assert (((operation == project_to_Q1) &&
                 (projection_matrix.m() == P) &&
                 (projection_matrix.n() == N) &&
                 (expansion_matrix.m() == N) &&
                 (expansion_matrix.n() == P))
                ||
                ((projection_matrix.m() == 0) &&
                 (projection_matrix.n() == 0)),
                ExcInternalError());

        switch (operation)
        {
          case none:
          {
            break;
          }
          case arithmetic_average:
          {
            double sum = 0;
            for (unsigned int i = 0; i < N; ++i)
              sum += values_out[i];

            const double average = sum / N;
            for (unsigned int i = 0; i < N; ++i)
              values_out[i] = average;

            break;
          }
          case harmonic_average:
          {
            // if one of the values is zero, the average is 0.0
            for (unsigned int i = 0; i < N; ++i)
              if (values_out[i] == 0.0)
              {
                for (unsigned int j = 0; j < N; ++j)
                  values_out[j] = 0.0;
                return;
              }

            double sum = 0;
            for (unsigned int i = 0; i < N; ++i)
              sum += 1. / values_out[i];

            const double average = 1. / (sum / N);
            for (unsigned int i = 0; i < N; ++i)
              values_out[i] = average;

            break;
          }
          case geometric_average:
          {
            double average = 1;
            for (unsigned int i = 0; i < N; ++i)
            {
              Assert (values_out[i] >= 0,
                      ExcMessage ("Computing the geometric average "
                                  "only makes sense for non-negative "
                                  "quantities."));
              average *= std::pow (values_out[i], 1. / N);
            }

            for (unsigned int i = 0; i < N; ++i)
              values_out[i] = average;

            break;
          }
          case pick_largest:
          {
            double max = -std::numeric_limits<double>::max();
            for (unsigned int i = 0; i < N; ++i)
              max = std::max(max, values_out[i]);

            for (unsigned int i = 0; i < N; ++i)
              values_out[i] = max;

            break;
          }
          case project_to_Q1:
          {
            // we will need the min/max values below, for use
            // after the projection operation
            double min =  std::numeric_limits<double>::max();
            double max = -std::numeric_limits<double>::max();
            for (unsigned int i = 0; i < N; ++i)
            {
              min = std::min(min, values_out[i]);
              max = std::max(max, values_out[i]);
            }

            // take the projection matrix and apply it to the values.
            Vector<double> x(N), z(P), y(N);
            for (unsigned int i = 0; i < N; ++i)
              y(i) = values_out[i];
            projection_matrix.vmult(z, y);

            // now that we have the Q1 values, restrict them to 
            // the min/max range of the original data
            for (unsigned int i = 0; i < P; ++i)
              z[i] = std::max (min, std::min (max, z[i]));

            // then expand back to the quadrature points
            expansion_matrix.vmult(x, z);
            for (unsigned int i = 0; i < N; ++i)
              values_out[i] = x(i);

            break;
          }
          case log_average:
          {
            double sum = 0;
            for (unsigned int i = 0; i < N; ++i)
            {
              if (values_out[i] == 0.0)
              {
                sum = -std::numeric_limits<double>::infinity();
                break;
              }
              Assert (values_out[i] > 0.0,
                      ExcMessage ("Computing the log average "
                                  "only makes sense for positive "
                                  "quantities."));
              sum += std::log10(values_out[i]);
            }

            const double log_value_average = std::pow(10., sum / N);
            for (unsigned int i = 0; i < N; ++i)
              values_out[i] = log_value_average;

            break;
          }
          default:
          {
            AssertThrow (false,
                         ExcMessage ("This averaging operation is not implemented."));
          }
        }

#ifdef DEBUG
#ifdef ASPECT_USE_FP_EXCEPTIONS
        // enable floating point exceptions again:
        feenableexcept(FE_DIVBYZERO|FE_INVALID);
#endif
#endif
      }


      template <int dim>
      void compute_projection_matrix (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                      const Quadrature<dim>     &quadrature_formula,
                                      const Mapping<dim>        &mapping,
                                      FullMatrix<double>        &projection_matrix,
                                      FullMatrix<double>        &expansion_matrix)
      {
        static FE_Q<dim> fe(1);
        FEValues<dim> fe_values (mapping, fe, quadrature_formula,
                                 update_values | update_JxW_values);

        const unsigned int P = fe.dofs_per_cell;
        const unsigned int N = quadrature_formula.size();

        FullMatrix<double> F(P, N);
        FullMatrix<double> M(P, P);

        projection_matrix.reinit(P, N);
        expansion_matrix.reinit(N, P);

        // reinitialize the fe_values object with the current cell. we get a
        // DoFHandler cell, but we are not going to use it with the
        // finite element associated with that DoFHandler, so cast it back
        // to just a tria iterator (all we need anyway is the geometry)
        fe_values.reinit (typename Triangulation<dim>::active_cell_iterator(cell));

        // compute the matrices F, M, E
        for (unsigned int i = 0; i < P; ++i)
          for (unsigned int q = 0; q < N; ++q)
            F(i, q) = fe_values.shape_value(i, q) * 
                      fe_values.JxW(q);

        for (unsigned int i = 0; i < P; ++i)
          for (unsigned int j = 0; j < P; ++j)
            for (unsigned int q = 0; q < N; ++q)
              M(i, j) += fe_values.shape_value(i, q) *
                         fe_values.shape_value(j, q) *
                         fe_values.JxW(q);

        for (unsigned int q = 0; q < N; ++q)
          for (unsigned int i = 0; i < P; ++i)
            expansion_matrix(q, i) = fe_values.shape_value(i, q);

        // replace M by M^{-1}
        M.gauss_jordan();

        // form M^{-1} F
        M.mmult(projection_matrix, F);
      }


      template <int dim>
      void average (const AveragingOperation operation,
                    const typename DoFHandler<dim>::active_cell_iterator &cell,
                    const Quadrature<dim> &quadrature_formula,
                    const Mapping<dim> &mapping,
                    const bool average_tangent_moduli,
                    MaterialModelOutputs<dim> &values_out)
      {
        FullMatrix<double> projection_matrix;
        FullMatrix<double> expansion_matrix;

        if (operation == project_to_Q1)
        {
          compute_projection_matrix (cell,
                                     quadrature_formula,
                                     mapping,
                                     projection_matrix,
                                     expansion_matrix);
        }

        if (average_tangent_moduli)
        {
          std::vector<double> tangent_modulus_components (values_out.tangent_moduli.size());
          for (unsigned int c = 0; c < SymmetricTensor<4,dim>::n_independent_components; ++c)
          {
            for (unsigned int i = 0; i < tangent_modulus_components.size(); ++i)
              tangent_modulus_components[i] = values_out.tangent_moduli[i].access_raw_entry(c);
            average_property (operation, projection_matrix, expansion_matrix,
                              tangent_modulus_components);
            for (unsigned int i = 0; i < tangent_modulus_components.size(); ++i)
              values_out.tangent_moduli[i].access_raw_entry(c) = tangent_modulus_components[i];
          }
        }

        average_property (operation, projection_matrix, expansion_matrix,
                          values_out.densities);
        average_property (operation, projection_matrix, expansion_matrix,
                          values_out.specific_heat);
        average_property (operation, projection_matrix, expansion_matrix,
                          values_out.thermal_expansivities);
        average_property (operation, projection_matrix, expansion_matrix,
                          values_out.thermal_conductivities);
      }
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace MaterialModel
  {
    namespace MaterialAveraging
    {
#define INSTANTIATE(dim) \
      template void \
      compute_projection_matrix (const typename DoFHandler<dim>::active_cell_iterator &, \
                                 const Quadrature<dim> &, \
                                 const Mapping<dim> &, \
                                 FullMatrix<double> &, \
                                 FullMatrix<double> &); \
      template void average (const AveragingOperation, \
                             const typename DoFHandler<dim>::active_cell_iterator &, \
                             const Quadrature<dim> &, \
                             const Mapping<dim> &, \
                             const bool, \
                             MaterialModelOutputs<dim> &);

      ELASPECT_INSTANTIATE(INSTANTIATE)
    }
  }
}
