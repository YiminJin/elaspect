#ifndef _elaspect_material_model_basic_properties_simple_h
#define _elaspect_material_model_basic_properties_simple_h

#include <elaspect/material_model/basic_properties/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace MaterialModel
  {
    namespace BasicProperties
    {
      using namespace dealii;

      template <int dim>
      class Simple : public Interface<dim>, public SimulatorAccess<dim>
      {
        public:
          void
          compute_basic_properties (const MaterialModel::MaterialModelInputs<dim> &in,
                                    const unsigned int i,
                                    MaterialModel::BasicPropertyOutputs<dim> &out) const override;

          FieldDependences::Dependence
          get_field_dependences (const MaterialProperties::Property requested_properties) const override;

          static
          void 
          declare_parameters (ParameterHandler &prm);

          virtual
          void
          parse_parameters (ParameterHandler &prm);

        private:
          double reference_temperature;
          std::vector<double> reference_densities;
          std::vector<double> specific_heat_capacities;
          std::vector<double> thermal_expansivities;
          std::vector<double> thermal_conductivities;
          std::vector<double> bulk_moduli;
          std::vector<double> shear_moduli;
      };
    }
  }
}

#endif
