#ifndef _elaspect_boundary_traction_zero_traction_h
#define _elaspect_boundary_traction_zero_traction_h

#include <elaspect/boundary_traction/interface.h>

namespace elaspect
{
  namespace BoundaryTraction
  {
    using namespace dealii;

    template <int dim>
    class ZeroTraction : public Interface<dim>
    {
      public:
        Tensor<1,dim>
        boundary_traction (const types::boundary_id boundary_indicator,
                           const Point<dim> &position,
                           const Tensor<1,dim> &normal_vector) const override;
    };
  }
}

#endif
