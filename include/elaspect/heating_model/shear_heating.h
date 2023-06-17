#ifndef _elaspect_heating_model_shear_heating_h
#define _elaspect_heating_model_shear_heating_h

#include <elaspect/heating_model/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace HeatingModel
  {
    using namespace dealii;

    template <int dim>
    class ShearHeating : public Interface<dim>, public SimulatorAccess<dim>
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
