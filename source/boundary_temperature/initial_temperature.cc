#include <elaspect/boundary_temperature/initial_temperature.h>
#include <elaspect/initial_temperature/interface.h>

namespace elaspect
{
  namespace BoundaryTemperature
  {
    template <int dim>
    double
    InitialTemperature<dim>::
    boundary_temperature(const types::boundary_id /*boundary_indicator*/,
                         const Point<dim> &position) const
    {
      return this->get_initial_temperature_manager().initial_temperature(position);
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace BoundaryTemperature
  {
    ELASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL(InitialTemperature,
                                                 "initial temperature",
                                                 "")
  }
}
