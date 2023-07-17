#ifndef _elaspect_postprocess_visualization_gravity_h
#define _elaspect_postprocess_visualization_gravity_h

#include <elaspect/postprocess/visualization.h>
#include <elaspect/simulator_access.h>

#include <deal.II/numerics/data_postprocessor.h>


namespace elaspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * A class derived from DataPostprocessorVector that outputs the gravity
       * as a vector field.
       */
      template <int dim>
      class Gravity
        : public DataPostprocessorVector<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          Gravity ();

          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double>> &computed_quantities) const override;
      };
    }
  }
}

#endif
