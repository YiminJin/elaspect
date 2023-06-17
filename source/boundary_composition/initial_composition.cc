#include <elaspect/boundary_composition/initial_composition.h>

namespace elaspect
{
  namespace BoundaryComposition
  {
    template <int dim>
    double
    InitialComposition<dim>::
    boundary_composition(const types::boundary_id /*boundary_indicator*/,
                         const Point<dim> &position,
                         const unsigned int compositional_field) const
    {
      return this->get_initial_composition_manager().initial_composition(position,
                                                                         compositional_field);
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace BoundaryComposition
  {
    ELASPECT_REGISTER_BOUNDARY_COMPOSITION_MODEL(InitialComposition,
                                             "initial composition",
                                             "")
  }
}
