#ifndef _elaspect_simulator_assemblers_interface_h
#define _elaspect_simulator_assemblers_interface_h

#include <elaspect/material_model/io_interface.h>
#include <elaspect/heating_model/interface.h>

#include <deal.II/fe/fe_values.h>

namespace elaspect
{
  using namespace dealii;

  namespace internal
  {
    namespace Assembly
    {
      namespace Scratch
      {
        template <int dim>
        struct ScratchBase
        {
          ScratchBase()
            : 
            cell(),
            face_number(numbers::invalid_unsigned_int)
          {}

          ScratchBase(const ScratchBase &scratch)
            : 
            cell(scratch.cell),
            face_number(scratch.face_number)
          {}

          virtual ~ScratchBase () = default;

          typename DoFHandler<dim>::active_cell_iterator cell;

          unsigned int face_number;
        };

        template <int dim>
        struct MechanicalSystem : public ScratchBase<dim>
        {
          MechanicalSystem (const FiniteElement<dim> &fe,
                       const Mapping<dim>       &mapping,
                       const Quadrature<dim>    &quadrature,
                       const Quadrature<dim-1>  &face_quadrature,
                       const UpdateFlags         update_flags,
                       const UpdateFlags         face_update_flags,
                       const unsigned int        u_dofs_per_cell,
                       const unsigned int        n_compositional_fields,
                       const MaterialModel::FieldDependences::Dependence field_dependences,
                       const MaterialModel::MaterialProperties::Property requested_properties);

          MechanicalSystem (const MechanicalSystem &scratch);

          FEValues<dim>     fe_values;
          FEFaceValues<dim> fe_face_values;

          std::vector<Tensor<1,dim> >          phi_u;
          std::vector<SymmetricTensor<2,dim> > symgrad_phi_u;

          std::vector<types::global_dof_index> local_dof_indices;

          MaterialModel::MaterialModelInputs<dim> material_model_inputs;
          MaterialModel::MaterialModelOutputs<dim> material_model_outputs;
        };

        template <int dim>
        struct ThermoSystem : public ScratchBase<dim>
        {
          ThermoSystem (const FiniteElement<dim> &fe,
                         const Mapping<dim>       &mapping,
                         const Quadrature<dim>    &quadrature,
                         const Quadrature<dim-1>  &face_quadrature,
                         const UpdateFlags         update_flags,
                         const UpdateFlags         face_update_flags,
                         const unsigned int        T_dofs_per_cell,
                         const unsigned int        n_compositional_fields,
                         const MaterialModel::FieldDependences::Dependence field_dependences,
                         const MaterialModel::MaterialProperties::Property requested_properties);

          ThermoSystem (const ThermoSystem &scratch);

          FEValues<dim>     fe_values;
          FEFaceValues<dim> fe_face_values;

          std::vector<double>           phi_T;
          std::vector<Tensor<1,dim>>    grad_phi_T;
          std::vector<double>           laplacian_phi_T;
          std::vector<double>           face_phi_T;

          std::vector<double>           old_temperature_values;
          std::vector<Tensor<1,dim>>    fluid_pressure_gradients;

          std::vector<Tensor<1,dim>>    displacement_increments;
          std::vector<Tensor<1,dim>>    mesh_displacement_increments;

          std::vector<types::global_dof_index> local_dof_indices;

          MaterialModel::MaterialModelInputs<dim>  material_model_inputs;
          MaterialModel::MaterialModelOutputs<dim> material_model_outputs;
          HeatingModel::HeatingModelOutputs        heating_model_outputs;
        };

        template <int dim>
        struct QPDSystem : public ScratchBase<dim>
        {
          QPDSystem (const FiniteElement<dim> &fe,
                     const Mapping<dim>       &mapping,
                     const Quadrature<dim>    &quadrature,
                     const Quadrature<dim-1>  &face_quadrature,
                     const UpdateFlags         update_flags,
                     const UpdateFlags         face_update_flags,
                     const unsigned int        field_dofs_per_cell,
                     const std::vector<typename Simulator<dim>::QPDField> &fields);

          QPDSystem (const QPDSystem &scratch);

          void reinit(const typename DoFHandler<dim>::active_cell_iterator &cell_ref);

          FEValues<dim> fe_values;

          std::unique_ptr<FEFaceValues<dim>>    fe_face_values;
          std::unique_ptr<FEFaceValues<dim>>    neighbor_fe_face_values;
          std::unique_ptr<FESubfaceValues<dim>> fe_subface_values;

          std::vector<types::global_dof_index>  local_dof_indices;
          std::vector<types::global_dof_index>  neighbor_dof_indices;

          std::vector<double>           phi_field;
          std::vector<Tensor<1,dim>>    grad_phi_field;
          std::vector<double>           face_phi_field;
          std::vector<double>           neighbor_face_phi_field;

          std::vector<Tensor<1,dim>>    displacement_increments;
          std::vector<Tensor<1,dim>>    mesh_displacement_increments;

          std::vector<Tensor<1,dim>>    face_displacement_increments;
          std::vector<Tensor<1,dim>>    face_mesh_displacement_increments;

          std::vector<std::vector<double>> old_field_values;

          const std::vector<typename Simulator<dim>::QPDField> qpd_fields;
        };
      }


      namespace CopyData
      {
        template <int dim>
        struct CopyDataBase
        {
          virtual ~CopyDataBase () = default;
        };

        template <int dim>
        struct MechanicalSystem : public CopyDataBase<dim>
        {
          MechanicalSystem (const unsigned int u_dofs_per_cell);

          MechanicalSystem (const MechanicalSystem &data);

          FullMatrix<double>    local_matrix;
          Vector<double>        local_rhs;

          std::vector<types::global_dof_index> local_dof_indices;
        };

        template <int dim>
        struct ThermoSystem : public CopyDataBase<dim>
        {
          ThermoSystem (const unsigned int T_dofs_per_cell);

          ThermoSystem (const ThermoSystem &data);

          FullMatrix<double> local_matrix;
          Vector<double>     local_rhs;

          std::vector<types::global_dof_index> local_dof_indices;
        };

        template <int dim>
        struct QPDSystem : public CopyDataBase<dim>
        {
          QPDSystem (const unsigned int field_dofs_per_cell,
                     const bool use_ALE_method,
                     const std::vector<typename Simulator<dim>::QPDField> &fields);

          FullMatrix<double> local_matrix;

          std::vector<FullMatrix<double>>   local_matrices_int_ext;
          std::vector<FullMatrix<double>>   local_matrices_ext_int;
          std::vector<FullMatrix<double>>   local_matrices_ext_ext;

          std::vector<bool>                 assembled_matrices;

          std::vector<Vector<double>>       local_rhs;

          std::vector<std::vector<types::global_dof_index>> local_dof_indices;
          std::vector<std::vector<types::global_dof_index>> neighbor_dof_indices;
        };
      }
    }
  }

  namespace Assemblers
  {
    template <int dim>
    class Interface
    {
      public:
        virtual ~Interface ();

        virtual
        void
        execute (internal::Assembly::Scratch::ScratchBase<dim> &scratch,
                 internal::Assembly::CopyData::CopyDataBase<dim> &data) const = 0;

        virtual
        void
        create_additional_material_model_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const;

        virtual
        MaterialModel::FieldDependences::Dependence get_field_dependences () const;

        virtual
        MaterialModel::MaterialProperties::Property get_needed_material_properties () const;
    };

    template <int dim>
    class Manager
    {
      public:
        void reset ();

        std::vector<std::unique_ptr<Assemblers::Interface<dim>>> thermo_system;

        std::vector<std::unique_ptr<Assemblers::Interface<dim>>> mechanical_system;

        std::vector<std::unique_ptr<Assemblers::Interface<dim>>> hydro_system;

        std::vector<std::unique_ptr<Assemblers::Interface<dim>>> qpd_system;

        std::vector<std::unique_ptr<Assemblers::Interface<dim>>> thermo_system_on_boundary_face;

        std::vector<std::unique_ptr<Assemblers::Interface<dim>>> mechanical_system_on_boundary_face;

        std::vector<std::unique_ptr<Assemblers::Interface<dim>>> qpd_system_on_boundary_face;

        std::vector<std::unique_ptr<Assemblers::Interface<dim>>> qpd_system_on_interior_face;
    };
  }
}

#endif
