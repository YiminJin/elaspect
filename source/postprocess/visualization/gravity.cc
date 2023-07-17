#include <elaspect/postprocess/visualization/gravity.h>
#include <elaspect/gravity_model/interface.h>

namespace elaspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      Gravity<dim>::Gravity()
        :
        DataPostprocessorVector<dim>("gravity", update_quadrature_points)
      {}


      template <int dim>
      void
      Gravity<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.evaluation_points.size();
        Assert (computed_quantities.size() == n_quadrature_points, ExcInternalError());
        Assert (computed_quantities[0].size() == dim, ExcInternalError());

        for (unsigned int q = 0; q < n_quadrature_points; ++q)
        {
          const Tensor<1,dim> g = this->get_gravity_model().gravity_vector(input_data.evaluation_points[q]);
          for (unsigned int k = 0; k < dim; ++k)
            computed_quantities[q](k) = g[k];
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
      ELASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(Gravity,
                                                    "gravity",
                                                    "A visualization output object that outputs the gravity vector."
                                                    "\n\n"
                                                    "Physical units: \\si {\\meter\\per\\second\\squared} .")
    }
  }
}
