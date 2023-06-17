#include <elaspect/boundary_traction/zero_traction.h>

namespace elaspect
{
  namespace BoundaryTraction
  {
    template <int dim>
    Tensor<1,dim>
    ZeroTraction<dim>::
    boundary_traction (const types::boundary_id,
                       const Point<dim> &,
                       const Tensor<1,dim> &) const
    {
      // return a zero tensor regardless of position
      return Tensor<1,dim>();
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace BoundaryTraction
  {
    ELASPECT_REGISTER_BOUNDARY_TRACTION_MODEL(ZeroTraction,
                                            "zero traction",
                                            "")
  }
}
