#ifndef _elaspect_simulator_assemblers_hydro_system_h
#define _elaspect_simulator_assemblers_hydro_system_h

#include <elaspect/simulator/assemblers/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace Assemblers
  {
    template <int dim>
    class HydroUnsaturatedFlow : public Interface<dim>, SimulatorAccess<dim>
    {
      public:
        void
        execute (internal::Assembly::Scratch::ScratchBase<dim> &scratch_base,
                 internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };
  }
}

#endif
