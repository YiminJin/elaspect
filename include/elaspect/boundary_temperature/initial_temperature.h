#ifndef _elaspect_boundary_temperature_initial_temperature_h
#define _elaspect_boundary_temperature_initial_temperature_h

#include <elaspect/boundary_temperature/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace BoundaryTemperature
  {
    template <int dim>
    class InitialTemperature : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        double boundary_temperature(const types::boundary_id boundary_indicator,
                                    const Point<dim> &position) const override;
    };
  }
}

#endif
