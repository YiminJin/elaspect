#ifndef _elaspect_time_stepping_inherited_time_step_h
#define _elaspect_time_stepping_inherited_time_step_h

#include <elaspect/time_stepping/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace TimeStepping
  {
    template <int dim>
    class InheritedTimeStep : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        double get_next_time_step_size() const override;
    };
  }
}

#endif
