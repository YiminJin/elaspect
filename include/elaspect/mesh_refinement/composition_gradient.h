#ifndef _elaspect_mesh_refinement_composition_gradient_h
#define _elaspect_mesh_refinement_composition_gradient_h

#include <elaspect/mesh_refinement/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace MeshRefinement
  {
    template <int dim>
    class CompositionGradient : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        void
        execute(Vector<float> &error_indicators) const override;

        static 
        void 
        declare_parameters(ParameterHandler &prm);

        void 
        parse_parameters(ParameterHandler &prm) override;

      private:
        std::vector<double> composition_scaling_factors;
    };
  }
}

#endif
