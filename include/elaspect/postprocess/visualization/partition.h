#ifndef _elaspect_postprocess_visualization_partition_h
#define _elaspect_postprocess_visualization_partition_h

#include <elaspect/postprocess/visualization.h>
#include <elaspect/simulator_access.h>

#include <deal.II/numerics/data_postprocessor.h>

namespace elaspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      class Partition
        : public DataPostprocessorScalar<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          Partition ();

          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double>> &computed_quantities) const override;
      };
    }
  }
}

#endif
