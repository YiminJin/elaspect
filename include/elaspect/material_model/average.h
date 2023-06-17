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
