#include <elaspect/postprocess/visualization/strain_rate.h>

namespace elaspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      StrainRate<dim>::StrainRate()
        :
        DataPostprocessorScalar<dim>("strain_rate", update_gradients)
      {}


      template <int dim>
      void
      StrainRate<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points, ExcInternalError());
        Assert (computed_quantities[0].size() == 1, ExcInternalError());
        Assert (input_data.solution_gradients[0].size() == this->introspection().n_components, ExcInternalError());

        for (unsigned int q = 0; q < n_quadrature_points; ++q)
        {
          Tensor<2,dim> grad_du;
          for (unsigned int d = 0; d < dim; ++d)
            grad_du[d] = input_data.solution_gradients[q][d];

          const SymmetricTensor<2,dim> depsilon = symmetrize(grad_du);
          computed_quantities[q](0) = std::sqrt(std::fabs(second_invariant(depsilon))) / this->get_timestep();
        }
      }
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ELASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(StrainRate,
                                                "strain rate",
                                                "")
    }
  }
}
