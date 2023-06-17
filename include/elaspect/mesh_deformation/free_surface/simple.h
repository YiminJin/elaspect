#ifndef _elaspect_mesh_deformation_free_surface_simple_h
#define _elaspect_mesh_deformation_free_surface_simple_h

#include <elaspect/mesh_deformation/free_surface/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace MeshDeformation
  {
    namespace FreeSurface
    {
      template <int dim>
      class Simple : public Interface<dim>, public SimulatorAccess<dim>
      {
        public:
          void
          make_boundary_constraints(const TrilinosWrappers::MPI::Vector &mesh_displacements,
                                    const DoFHandler<dim>               &mesh_deformation_dof_handler,
                                    AffineConstraints<double>           &mesh_deformation_constraints,
                                    const std::set<types::boundary_id>  &boundary_id) const override;
      };
    }
  }
}

#endif
