#include <elaspect/initial_topography/zero_topography.h>

namespace elaspect
{
  namespace InitialTopography
  {
    template <int dim>
    double
    ZeroTopography<dim>::
    value (const Point<dim-1> &/*p*/) const
    {
      // return a zero value regardless of position
      return 0.0;
    }


    template <int dim>
    double
    ZeroTopography<dim>::
    max_topography () const
    {
      return 0;
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace InitialTopography
  {
    ELASPECT_REGISTER_INITIAL_TOPOGRAPHY_MODEL(ZeroTopography,
                                           "zero topography",
                                           "Implementation of a model in which the initial topography "
                                           "is zero. ")
  }
}
