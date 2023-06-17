#ifndef _elaspect_material_model_viscosity_simple_h
#define _elaspect_material_model_viscosity_simple_h

#include <elaspect/material_model/viscosity/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace MaterialModel
  {
    namespace Viscosity
    {
      using namespace dealii;

      template <int dim>
      class Simple : public Interface<dim>, public SimulatorAccess<dim>
      {
        public:
          double
          compute_viscosity (const MaterialModel::MaterialModelInputs<dim> &in,
                             const unsigned int i,
                             const unsigned int field) const override;

          FieldDependences::Dependence
          get_field_dependences () const override;

          static
          void
          declare_parameters (ParameterHandler &prm);

          virtual
          void
          parse_parameters (ParameterHandler &prm);

        private:
          double reference_T;
          std::vector<double> eta;
          std::vector<double> beta;
          double maximum_thermal_prefactor;
          double minimum_thermal_prefactor;
      };
    }
  }
}

#endif
