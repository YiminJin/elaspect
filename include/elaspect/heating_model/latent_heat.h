#ifndef _elaspect_heating_model_latent_heat_h
#define _elaspect_heating_model_latent_heat_h

#include <elaspect/heating_model/interface.h>

namespace elaspect
{
  namespace HeatingModel
  {
    using namespace dealii;

    template <int dim>
    class LatentHeat : public Interface<dim>
    {
      public:
        /**
         * Compute the heating model outputs for this class.
         */
        void
        evaluate(const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                 const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                 HeatingModel::HeatingModelOutputs &heating_model_outputs) const override;
    };
  }
}

#endif
