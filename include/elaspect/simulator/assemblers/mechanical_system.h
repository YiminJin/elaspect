#ifndef _elaspect_simulator_assemblers_mechanical_system_h
#define _elaspect_simulator_assemblers_mechanical_system_h

#include <elaspect/simulator/assemblers/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace Assemblers
  {
    template <int dim>
    class MechanicalQuasiStaticTerms : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        void
        execute (internal::Assembly::Scratch::ScratchBase<dim> &scratch_base,
                 internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;

        MaterialModel::FieldDependences::Dependence
        get_field_dependences () const override;

        MaterialModel::MaterialProperties::Property
        get_needed_material_properties () const override;
    };

    template <int dim>
    class MechanicalBoundaryTraction : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        void
        execute (internal::Assembly::Scratch::ScratchBase<dim> &scratch_base,
                 internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };
  }
}

#endif
