#ifndef _elaspect_postprocess_visualization_strain_rate_h
#define _elaspect_postprocess_visualization_strain_rate_h

#include <elaspect/postprocess/visualization.h>

#include <deal.II/numerics/data_postprocessor.h>

namespace elaspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      class StrainRate
        : public DataPostprocessorScalar<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          /**
           * Constructor.
           */
          StrainRate();

          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double> > &computed_quantities) const override;
      };
    }
  }
}

#endif
