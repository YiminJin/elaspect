#include <elaspect/simulator/assemblers/interface.h>
#include <elaspect/simulator.h>

#include <deal.II/base/signaling_nan.h>

namespace elaspect
{
  namespace internal
  {
    namespace Assembly
    {
      namespace Scratch
      {
        template <int dim>
        MechanicalSystem<dim>::
        MechanicalSystem(const FiniteElement<dim> &fe,
                         const Mapping<dim>       &mapping,
                         const Quadrature<dim>    &quadrature,
                         const Quadrature<dim-1>  &face_quadrature,
                         const UpdateFlags         update_flags,
                         const UpdateFlags         face_update_flags,
                         const unsigned int        u_dofs_per_cell,
                         const unsigned int        n_compositional_fields,
                         const MaterialModel::FieldDependences::Dependence field_dependences,
                         const MaterialModel::MaterialProperties::Property requested_properties)
          :
          ScratchBase<dim>(),

          fe_values (mapping, fe, quadrature, update_flags),
          fe_face_values (mapping, fe, face_quadrature, face_update_flags),

          phi_u (u_dofs_per_cell, numbers::signaling_nan<Tensor<1,dim> >()),
          symgrad_phi_u (u_dofs_per_cell, numbers::signaling_nan<SymmetricTensor<2,dim> >()),

          local_dof_indices (fe.dofs_per_cell),

          material_model_inputs (quadrature.size(), n_compositional_fields,
                                 field_dependences, requested_properties),
          material_model_outputs (quadrature.size(), n_compositional_fields)
        {}


        template <int dim>
        MechanicalSystem<dim>::MechanicalSystem (const MechanicalSystem &scratch)
          :
          ScratchBase<dim>(scratch),

          fe_values (scratch.fe_values.get_mapping(),
                     scratch.fe_values.get_fe(),
                     scratch.fe_values.get_quadrature(),
                     scratch.fe_values.get_update_flags()),
          fe_face_values (scratch.fe_face_values.get_mapping(),
                          scratch.fe_face_values.get_fe(),
                          scratch.fe_face_values.get_quadrature(),
                          scratch.fe_face_values.get_update_flags()),

          phi_u (scratch.phi_u),
          symgrad_phi_u (scratch.symgrad_phi_u),
          local_dof_indices (scratch.local_dof_indices),
          material_model_inputs (scratch.material_model_inputs),
          material_model_outputs (scratch.material_model_outputs)
        {}


        template <int dim>
        ThermoSystem<dim>::
        ThermoSystem (const FiniteElement<dim> &fe,
                      const Mapping<dim>       &mapping,
                      const Quadrature<dim>    &quadrature,
                      const Quadrature<dim-1>  &face_quadrature,
                      const UpdateFlags         update_flags,
                      const UpdateFlags         face_update_flags,
                      const unsigned int        T_dofs_per_cell,
                      const unsigned int        n_compositional_fields,
                      const MaterialModel::FieldDependences::Dependence field_dependences,
                      const MaterialModel::MaterialProperties::Property requested_properties)
          :
          ScratchBase<dim>(),

          fe_values(mapping, fe, quadrature, update_flags),
          fe_face_values(face_quadrature.size() > 0
                         ?
                         std::make_unique<FEFaceValues<dim>>(mapping, fe, face_quadrature, face_update_flags)
                         :
                         nullptr),
          neighbor_fe_face_values(face_quadrature.size() > 0
                                  ?
                                  std::make_unique<FEFaceValues<dim>>(mapping, fe, face_quadrature, face_update_flags)
                                  :
                                  nullptr),
          fe_subface_values(face_quadrature.size() > 0
                            ?
                            std::make_unique<FESubfaceValues<dim>>(mapping, fe, face_quadrature, face_update_flags)
                            :
                            nullptr),

          local_dof_indices(fe.dofs_per_cell),
          neighbor_dof_indices(fe.dofs_per_cell),

          phi_T(T_dofs_per_cell, numbers::signaling_nan<double>()),
          grad_phi_T(T_dofs_per_cell, numbers::signaling_nan<Tensor<1,dim>>()),
          face_phi_T(face_quadrature.size() > 0 ? T_dofs_per_cell : 0,
                     numbers::signaling_nan<double>()),
          face_grad_phi_T(face_quadrature.size() > 0 ? T_dofs_per_cell : 0,
                          numbers::signaling_nan<Tensor<1,dim>>()),
          neighbor_face_phi_T(face_quadrature.size() > 0 ? T_dofs_per_cell : 0,
                              numbers::signaling_nan<double>()),
          neighbor_face_grad_phi_T(face_quadrature.size() > 0 ? T_dofs_per_cell : 0,
                                   numbers::signaling_nan<Tensor<1,dim>>()),

          old_temperature_values(quadrature.size(), numbers::signaling_nan<double>()),
          fluid_pressure_gradients(quadrature.size(), numbers::signaling_nan<Tensor<1,dim>>()),

          displacement_increments(quadrature.size(), numbers::signaling_nan<Tensor<1,dim>>()),
          mesh_displacement_increments(quadrature.size(), numbers::signaling_nan<Tensor<1,dim>>()),
          face_displacement_increments(face_quadrature.size(), numbers::signaling_nan<Tensor<1,dim>>()),
          face_mesh_displacement_increments(face_quadrature.size(), numbers::signaling_nan<Tensor<1,dim>>()),

          material_model_inputs(quadrature.size(), n_compositional_fields,
                                field_dependences, requested_properties),
          material_model_outputs(quadrature.size(), n_compositional_fields),

          face_material_model_inputs(face_quadrature.size(), n_compositional_fields,
                                     field_dependences, requested_properties),
          face_material_model_outputs(face_quadrature.size(), n_compositional_fields),

          neighbor_face_material_model_inputs(face_quadrature.size(), n_compositional_fields,
                                              field_dependences, requested_properties),
          neighbor_face_material_model_outputs(face_quadrature.size(), n_compositional_fields),

          heating_model_outputs(quadrature.size()),
          face_heating_model_outputs(face_quadrature.size()),
          neighbor_face_heating_model_outputs(face_quadrature.size())
        {}


        template <int dim>
        ThermoSystem<dim>::ThermoSystem (const ThermoSystem &scratch)
          :
          ScratchBase<dim>(scratch),

          fe_values(scratch.fe_values.get_mapping(),
                    scratch.fe_values.get_fe(),
                    scratch.fe_values.get_quadrature(),
                    scratch.fe_values.get_update_flags()),
          fe_face_values(scratch.fe_face_values.get()
                         ?
                         std::make_unique<FEFaceValues<dim>>(scratch.fe_face_values->get_mapping(),
                                                             scratch.fe_face_values->get_fe(),
                                                             scratch.fe_face_values->get_quadrature(),
                                                             scratch.fe_face_values->get_update_flags())
                         :
                         nullptr),
          neighbor_fe_face_values(scratch.neighbor_fe_face_values.get()
                                  ?
                                  std::make_unique<FEFaceValues<dim>>(scratch.neighbor_fe_face_values->get_mapping(),
                                                                      scratch.neighbor_fe_face_values->get_fe(),
                                                                      scratch.neighbor_fe_face_values->get_quadrature(),
                                                                      scratch.neighbor_fe_face_values->get_update_flags())
                                  :
                                  nullptr),
          fe_subface_values(scratch.fe_subface_values.get()
                            ?
                            std::make_unique<FESubfaceValues<dim>>(scratch.fe_subface_values->get_mapping(),
                                                                   scratch.fe_subface_values->get_fe(),
                                                                   scratch.fe_subface_values->get_quadrature(),
                                                                   scratch.fe_subface_values->get_update_flags())
                            :
                            nullptr),

          local_dof_indices(scratch.local_dof_indices),
          neighbor_dof_indices(scratch.neighbor_dof_indices),

          phi_T(scratch.phi_T),
          grad_phi_T(scratch.grad_phi_T),
          face_phi_T(scratch.face_phi_T),
          face_grad_phi_T(scratch.face_grad_phi_T),
          neighbor_face_phi_T(scratch.neighbor_face_phi_T),
          neighbor_face_grad_phi_T(scratch.neighbor_face_grad_phi_T),

          old_temperature_values(scratch.old_temperature_values),
          fluid_pressure_gradients(scratch.fluid_pressure_gradients),

          displacement_increments(scratch.displacement_increments),
          mesh_displacement_increments(scratch.mesh_displacement_increments),
          face_displacement_increments(scratch.face_displacement_increments),
          face_mesh_displacement_increments(scratch.face_mesh_displacement_increments),

          material_model_inputs(scratch.material_model_inputs),
          material_model_outputs(scratch.material_model_outputs),
          face_material_model_inputs(scratch.face_material_model_inputs),
          face_material_model_outputs(scratch.face_material_model_outputs),
          neighbor_face_material_model_inputs(scratch.neighbor_face_material_model_inputs),
          neighbor_face_material_model_outputs(scratch.neighbor_face_material_model_outputs),

          heating_model_outputs(scratch.heating_model_outputs),
          face_heating_model_outputs(scratch.face_heating_model_outputs),
          neighbor_face_heating_model_outputs(scratch.neighbor_face_heating_model_outputs)
        {}


        template <int dim>
        QPDSystem<dim>::QPDSystem(const FiniteElement<dim> &fe,
                                  const Mapping<dim>       &mapping,
                                  const Quadrature<dim>    &quadrature,
                                  const Quadrature<dim-1>  &face_quadrature,
                                  const UpdateFlags         update_flags,
                                  const UpdateFlags         face_update_flags,
                                  const unsigned int        field_dofs_per_cell,
                                  const std::vector<typename Simulator<dim>::QPDField> &fields)
          :
          ScratchBase<dim>(),

          fe_values(mapping, fe, quadrature, update_flags),
          fe_face_values(face_quadrature.size() > 0
                         ?
                         std::make_unique<FEFaceValues<dim>>(mapping, fe, face_quadrature, face_update_flags)
                         :
                         nullptr),
          neighbor_fe_face_values(face_quadrature.size() > 0
                                  ?
                                  std::make_unique<FEFaceValues<dim>>(mapping, fe, face_quadrature, face_update_flags)
                                  :
                                  nullptr),
          fe_subface_values(face_quadrature.size() > 0
                            ?
                            std::make_unique<FESubfaceValues<dim>>(mapping, fe, face_quadrature, face_update_flags)
                            :
                            nullptr),

          local_dof_indices(fe.dofs_per_cell),
          neighbor_dof_indices(fe.dofs_per_cell),

          phi_field(field_dofs_per_cell, numbers::signaling_nan<double>()),
          grad_phi_field(field_dofs_per_cell, numbers::signaling_nan<Tensor<1,dim>>()),
          face_phi_field(face_quadrature.size() > 0 ? field_dofs_per_cell : 0, 
                         numbers::signaling_nan<double>()),
          neighbor_face_phi_field(face_quadrature.size() > 0 ? field_dofs_per_cell : 0, 
                                  numbers::signaling_nan<double>()),

          displacement_increments(quadrature.size(), numbers::signaling_nan<Tensor<1,dim>>()),
          mesh_displacement_increments(quadrature.size(), numbers::signaling_nan<Tensor<1,dim>>()),
          face_displacement_increments(face_quadrature.size(), numbers::signaling_nan<Tensor<1,dim>>()),
          face_mesh_displacement_increments(face_quadrature.size(), numbers::signaling_nan<Tensor<1,dim>>()),
          old_field_values(fields.size(), 
                           std::vector<double>(quadrature.size(), numbers::signaling_nan<double>())),
          qpd_fields(fields)
        {}


        template <int dim>
        QPDSystem<dim>::
        QPDSystem(const QPDSystem &scratch)
          :
          ScratchBase<dim>(scratch),

          fe_values(scratch.fe_values.get_mapping(),
                    scratch.fe_values.get_fe(),
                    scratch.fe_values.get_quadrature(),
                    scratch.fe_values.get_update_flags()),
          fe_face_values(scratch.fe_face_values.get()
                         ?
                         std::make_unique<FEFaceValues<dim>>(scratch.fe_face_values->get_mapping(),
                                                             scratch.fe_face_values->get_fe(),
                                                             scratch.fe_face_values->get_quadrature(),
                                                             scratch.fe_face_values->get_update_flags())
                         :
                         nullptr),
          neighbor_fe_face_values(scratch.neighbor_fe_face_values.get()
                                  ?
                                  std::make_unique<FEFaceValues<dim>>(scratch.neighbor_fe_face_values->get_mapping(),
                                                                      scratch.neighbor_fe_face_values->get_fe(),
                                                                      scratch.neighbor_fe_face_values->get_quadrature(),
                                                                      scratch.neighbor_fe_face_values->get_update_flags())
                                  :
                                  nullptr),
          fe_subface_values(scratch.fe_subface_values.get()
                            ?
                            std::make_unique<FESubfaceValues<dim>>(scratch.fe_subface_values->get_mapping(),
                                                                   scratch.fe_subface_values->get_fe(),
                                                                   scratch.fe_subface_values->get_quadrature(),
                                                                   scratch.fe_subface_values->get_update_flags())
                            :
                            nullptr),

          local_dof_indices(scratch.local_dof_indices),
          neighbor_dof_indices(scratch.neighbor_dof_indices),

          phi_field(scratch.phi_field),
          grad_phi_field(scratch.grad_phi_field),
          face_phi_field(scratch.face_phi_field),
          neighbor_face_phi_field(scratch.neighbor_face_phi_field),
          
          displacement_increments(scratch.displacement_increments),
          mesh_displacement_increments(scratch.mesh_displacement_increments),
          face_displacement_increments(scratch.face_displacement_increments),
          face_mesh_displacement_increments(scratch.face_mesh_displacement_increments),
          old_field_values(scratch.old_field_values),

          qpd_fields(scratch.qpd_fields)
        {}
      }


      namespace CopyData
      {
        template <int dim>
        MechanicalSystem<dim>::MechanicalSystem (const unsigned int u_dofs_per_cell)
          :
          local_matrix (u_dofs_per_cell, u_dofs_per_cell),
          local_rhs (u_dofs_per_cell),
          local_dof_indices (u_dofs_per_cell)
        {}


        template <int dim>
        MechanicalSystem<dim>::MechanicalSystem (const MechanicalSystem<dim> &data)
          :
          local_matrix (data.local_matrix),
          local_rhs (data.local_rhs),
          local_dof_indices (data.local_dof_indices)
        {}


        template <int dim>
        ThermoSystem<dim>::ThermoSystem (const unsigned int T_dofs_per_cell,
                                         const bool is_discontinuous)
          :
          local_matrix (T_dofs_per_cell, T_dofs_per_cell),

          local_matrices_int_ext((is_discontinuous
                                  ?
                                  GeometryInfo<dim>::max_children_per_face * GeometryInfo<dim>::faces_per_cell
                                  :
                                  0),
                                 FullMatrix<double>(T_dofs_per_cell, T_dofs_per_cell)),
          local_matrices_ext_int((is_discontinuous
                                  ?
                                  GeometryInfo<dim>::max_children_per_face * GeometryInfo<dim>::faces_per_cell
                                  :
                                  0),
                                 FullMatrix<double>(T_dofs_per_cell, T_dofs_per_cell)),
          local_matrices_ext_ext((is_discontinuous
                                  ?
                                  GeometryInfo<dim>::max_children_per_face * GeometryInfo<dim>::faces_per_cell
                                  :
                                  0),
                                 FullMatrix<double>(T_dofs_per_cell, T_dofs_per_cell)),

          assembled_matrices((is_discontinuous
                              ?
                              GeometryInfo<dim>::max_children_per_face * GeometryInfo<dim>::faces_per_cell
                              :
                              0), false),

          local_rhs(T_dofs_per_cell),
          local_dof_indices(T_dofs_per_cell),
          neighbor_dof_indices((is_discontinuous
                                ?
                                GeometryInfo<dim>::max_children_per_face * GeometryInfo<dim>::faces_per_cell
                                :
                                0), std::vector<types::global_dof_index>(T_dofs_per_cell))
        {}


        template <int dim>
        QPDSystem<dim>::QPDSystem(const unsigned int field_dofs_per_cell,
                                  const bool is_discontinuous,
                                  const std::vector<typename Simulator<dim>::QPDField> &fields)
          :
          local_matrix(field_dofs_per_cell, field_dofs_per_cell),

          local_matrices_int_ext((is_discontinuous
                                  ?
                                  GeometryInfo<dim>::max_children_per_face * GeometryInfo<dim>::faces_per_cell
                                  :
                                  0),
                                 FullMatrix<double>(field_dofs_per_cell, field_dofs_per_cell)),
          local_matrices_ext_int((is_discontinuous
                                  ?
                                  GeometryInfo<dim>::max_children_per_face * GeometryInfo<dim>::faces_per_cell
                                  :
                                  0),
                                 FullMatrix<double>(field_dofs_per_cell, field_dofs_per_cell)),
          local_matrices_ext_ext((is_discontinuous
                                  ?
                                  GeometryInfo<dim>::max_children_per_face * GeometryInfo<dim>::faces_per_cell
                                  :
                                  0),
                                 FullMatrix<double>(field_dofs_per_cell, field_dofs_per_cell)),

          assembled_matrices((is_discontinuous
                              ?
                              GeometryInfo<dim>::max_children_per_face * GeometryInfo<dim>::faces_per_cell
                              :
                              0), false),

          local_rhs(fields.size(), Vector<double>(field_dofs_per_cell)),
          local_dof_indices(fields.size(), std::vector<types::global_dof_index>(field_dofs_per_cell)),

          neighbor_dof_indices((is_discontinuous
                                ?
                                GeometryInfo<dim>::max_children_per_face * GeometryInfo<dim>::faces_per_cell
                                :
                                0), std::vector<types::global_dof_index>(field_dofs_per_cell))
        {}
      }
    }
  }


  namespace Assemblers
  {
    template <int dim>
    Interface<dim>::~Interface ()
    {}


    template <int dim>
    void
    Interface<dim>::
    create_additional_material_model_outputs (MaterialModel::MaterialModelOutputs<dim> &) const
    {}


    template <int dim>
    MaterialModel::FieldDependences::Dependence
    Interface<dim>::get_field_dependences () const
    {
      return MaterialModel::FieldDependences::none;
    }


    template <int dim>
    MaterialModel::MaterialProperties::Property
    Interface<dim>::get_needed_material_properties () const
    {
      return MaterialModel::MaterialProperties::none;
    }


    template <int dim>
    void Manager<dim>::reset ()
    {
      mechanical_system.clear();
      mechanical_system_on_boundary_face.clear();
      qpd_system.clear();
      qpd_system_on_boundary_face.clear();
      qpd_system_on_interior_face.clear();
      thermo_system.clear();
      thermo_system_on_boundary_face.clear();
      thermo_system_on_interior_face.clear();
      hydro_system.clear();
    }
  }
}


// explicit instantiations
namespace elaspect
{
#define INSTANTIATE(dim) \
  namespace internal { \
    namespace Assembly { \
      namespace Scratch { \
        template struct MechanicalSystem<dim>; \
        template struct ThermoSystem<dim>; \
        template struct QPDSystem<dim>; \
      } \
      namespace CopyData { \
        template struct MechanicalSystem<dim>; \
        template struct ThermoSystem<dim>; \
        template struct QPDSystem<dim>; \
      } \
    } \
  } \
  namespace Assemblers { \
    template class Interface<dim>; \
    template class Manager<dim>; \
  }

  ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
