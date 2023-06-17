#include <elaspect/mesh_refinement/velocity.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/numerics/error_estimator.h>

namespace elaspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    Velocity<dim>::execute (Vector<float> &indicators) const
    {
      const Introspection<dim> &introspection = this->introspection();

      TrilinosWrappers::MPI::BlockVector 
      disp_incr (introspection.index_sets.system_partitioning,
                 introspection.index_sets.system_relevant_partitioning,
                 this->get_mpi_communicator()),
      dist_disp_incr (introspection.index_sets.system_partitioning,
                      this->get_mpi_communicator()),
      dist_disp_old (introspection.index_sets.system_partitioning,
                     this->get_mpi_communicator());

      const unsigned int block_idx = introspection.block_indices.displacement;
      dist_disp_incr.block(block_idx) = this->get_solution().block(block_idx);
      dist_disp_old.block(block_idx)  = this->get_old_solution().block(block_idx);
      dist_disp_incr.block(block_idx).add (-1, dist_disp_old.block(block_idx));
      disp_incr.block(block_idx) = dist_disp_old.block(block_idx);

      indicators = 0;

      const QGauss<dim-1> quadrature (this->get_parameters().displacement_degree + 1);
      KellyErrorEstimator<dim>::estimate (this->get_mapping(),
                                          this->get_dof_handler(),
                                          quadrature,
                                          std::map<types::boundary_id, const Function<dim>*>(),
                                          disp_incr,
                                          indicators,
                                          introspection.component_masks.displacement,
                                          nullptr,
                                          0,
                                          this->get_triangulation().locally_owned_subdomain());
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace MeshRefinement
  {
    ELASPECT_REGISTER_MESH_REFINEMENT_CRITERION(Velocity,
                                            "velocity",
                                            "")
  }
}
