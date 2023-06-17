#ifndef _elaspect_material_model_viscosity_diffusion_dislocation_h
#define _elaspect_material_model_viscosity_diffusion_dislocation_h

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
      class DiffusionDislocation : public Interface<dim>, public SimulatorAccess<dim>
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
          std::vector<double> A_diff;
          std::vector<double> m_diff;
          std::vector<double> E_diff;
          std::vector<double> V_diff;

          std::vector<double> A_disl;
          std::vector<double> n_disl;
          std::vector<double> E_disl;
          std::vector<double> V_disl;

          double d;
          double eta_min;
          double eta_max;
      };
    }
  }
}

#endif
