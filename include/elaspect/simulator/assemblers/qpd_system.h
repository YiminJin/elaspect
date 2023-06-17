#ifndef _elaspect_simulator_assemblers_qpd_system_h
#define _elaspect_simulator_assemblers_qpd_system_h

#include <elaspect/simulator/assemblers/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace Assemblers
  {
    template <int dim>
    class QPDProjection : public Assemblers::Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };

    template <int dim>
    class QPDAdvection : public Assemblers::Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };

    template <int dim>
    class QPDAdvectionBoundaryFace : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };

    template <int dim>
    class QPDAdvectionInteriorFace : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };
  }
}

#endif
