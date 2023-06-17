#ifndef _elaspect_mesh_deformation_mesh_smoothing_laplacian_h
#define _elaspect_mesh_deformation_mesh_smoothing_laplacian_h

#include <elaspect/mesh_deformation/mesh_smoothing/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace MeshDeformation
  {
    namespace MeshSmoothing
    {
      template <int dim>
      class Laplacian : public Interface<dim>, public SimulatorAccess<dim>
      {
        public:
          void
          compute_mesh_displacements(const DoFHandler<dim>           &dof_handler,
                                     const AffineConstraints<double> &constraints,
                                     TrilinosWrappers::MPI::Vector   &displacements) const override;
      };
    }
  }
}

#endif
