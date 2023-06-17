#ifndef _elaspect_initial_topography_zero_topography_h
#define _elaspect_initial_topography_zero_topography_h

#include <elaspect/initial_topography/interface.h>

namespace elaspect
{
  namespace InitialTopography
  {
    template <int dim>
    class ZeroTopography : public Interface<dim>
    {
      public:
        double
        value (const Point<dim-1> &p) const override;

        double max_topography () const override;
    };
  }
}

#endif
