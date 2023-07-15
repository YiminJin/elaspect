#include <elaspect/time_stepping/inherited_time_step.h>

namespace elaspect
{
  namespace TimeStepping
  {
    template <int dim>
    double
    InheritedTimeStep<dim>::get_next_time_step_size() const
    {
      if (this->get_timestep_number() == 0)
        return std::numeric_limits<double>::max();
      else
        return this->get_timestep();
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace TimeStepping
  {
    ELASPECT_REGISTER_TIME_STEPPING_MODEL(InheritedTimeStep,
                                          "inherited time step",
                                          "")
  }
}
