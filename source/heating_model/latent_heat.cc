#include <elaspect/heating_model/latent_heat.h>

namespace elaspect
{
  namespace HeatingModel
  {
    template <int dim>
    void
    LatentHeat<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
             const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
             HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      const unsigned int n_evaluation_points = material_model_inputs.n_evaluation_points();
      AssertDimension(heating_model_outputs.heating_source_terms.size(), n_evaluation_points);

      for (unsigned int q = 0; q < n_evaluation_points; ++q)
      {
        heating_model_outputs.heating_source_terms[q] = 
          material_model_inputs.temperature[q] * material_model_outputs.densities[q];

        heating_model_outputs.lhs_latent_heat_terms[q] =
          - material_model_inputs.temperature[q] * material_model_outputs.densities[q]
          * material_model_outputs.entropy_derivative_temperature[q];
      }
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace HeatingModel
  {
    ELASPECT_REGISTER_HEATING_MODEL(LatentHeat,
                                "latent heat",
                                "Implementation of a standard model for latent heat.")
  }
}
