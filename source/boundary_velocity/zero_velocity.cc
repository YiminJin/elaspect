#include <elaspect/boundary_velocity/zero_velocity.h>


namespace elaspect
{
  namespace BoundaryVelocity
  {
    template <int dim>
    Tensor<1,dim>
    ZeroVelocity<dim>::
    boundary_velocity (const types::boundary_id,
                       const Point<dim> &) const
    {
      // return a zero tensor regardless of position
      return Tensor<1,dim>();
    }
  }
}

// explicit instantiations
namespace elaspect
{
  namespace BoundaryVelocity
  {
    ELASPECT_REGISTER_BOUNDARY_VELOCITY_MODEL(ZeroVelocity,
                                          "zero velocity",
                                          "Implementation of a model in which the boundary "
                                          "velocity is zero. This is commonly referred to as "
                                          "a ``stick boundary condition'', indicating that "
                                          "the material ``sticks'' to the material on the "
                                          "other side of the boundary.")
  }
}
