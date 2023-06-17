#ifndef _elaspect_boundary_velocity_zero_velocity_h
#define _elaspect_boundary_velocity_zero_velocity_h

#include <elaspect/boundary_velocity/interface.h>

namespace elaspect
{
  namespace BoundaryVelocity
  {
    using namespace dealii;

    template <int dim>
    class ZeroVelocity : public Interface<dim>
    {
      public:
        /**
         * Return the boundary velocity as a function of position. For the
         * current class, this function obviously simply returns a zero
         * tensor.
         */
        Tensor<1,dim>
        boundary_velocity (const types::boundary_id boundary_indicator,
                           const Point<dim> &position) const override;

        // avoid -Woverloaded-virtual warning until the deprecated function
        // is removed from the interface:
        using Interface<dim>::boundary_velocity;
    };
  }
}


#endif
