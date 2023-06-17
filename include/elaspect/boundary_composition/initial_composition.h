#ifndef _elaspect_boundary_composition_initial_composition_h
#define _elaspect_boundary_composition_initial_composition_h

#include <elaspect/initial_composition/interface.h>
#include <elaspect/boundary_composition/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace BoundaryComposition
  {
    template <int dim>
    class InitialComposition : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        double boundary_composition (const types::boundary_id boundary_indicator,
                                     const Point<dim> &position,
                                     const unsigned int compositional_field) const override;
    };
  }
}

#endif
