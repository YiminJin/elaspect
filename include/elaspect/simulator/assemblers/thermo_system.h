#ifndef _elaspect_simulator_assemblers_thermo_system_h
#define _elaspect_simulator_assemblers_thermo_system_h

#include <elaspect/simulator/assemblers/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace Assemblers
  {
    template <int dim>
    class ThermoConvectionDiffusion : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;

        MaterialModel::MaterialProperties::Property
        get_needed_material_properties() const override;
    };

    template <int dim>
    class ThermoBoundaryFlux : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };

    template <int dim>
    class ThermoSystemBoundaryFace : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;

        MaterialModel::MaterialProperties::Property
        get_needed_material_properties() const override;
    };

    template <int dim>
    class ThermoSystemInteriorFace : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;

        MaterialModel::MaterialProperties::Property
        get_needed_material_properties() const override;
    };
  }
}

#endif
