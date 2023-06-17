#ifndef _elaspect_mesh_refinement_velocity_h
#define _elaspect_mesh_refinement_velocity_h

#include <elaspect/mesh_refinement/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace MeshRefinement
  {
    template <int dim>
    class Velocity : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        void execute (Vector<float> &error_indicators) const override;
    };
  }
}

#endif
