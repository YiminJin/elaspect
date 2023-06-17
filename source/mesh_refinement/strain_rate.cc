#include <elaspect/mesh_refinement/strain_rate.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace elaspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    StrainRate<dim>::execute (Vector<float> &indicators) const
    {
      indicators = 0;

      const QMidpoint<dim> quadrature;

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature,
                               update_quadrature_points | update_values | update_gradients);

      std::vector<SymmetricTensor<2,dim> > strain_increments (quadrature.size());

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
        {
          const unsigned int idx = cell->active_cell_index();
          fe_values.reinit(cell);

          fe_values[this->introspection().extractors.displacement].get_function_symmetric_gradients (
            this->get_solution(), strain_increments);

          indicators(idx) = strain_increments[0].norm();
        }
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace MeshRefinement
  {
    ELASPECT_REGISTER_MESH_REFINEMENT_CRITERION(StrainRate,
                                            "strain rate",
                                            "A mesh refinement criterion that computes the "
                                            "refinement indicators equal to the strain rate "
                                            "norm computed at the center of the elements.")
  }
}
