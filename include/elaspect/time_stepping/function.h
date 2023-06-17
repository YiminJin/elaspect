#ifndef _elaspect_time_stepping_function_h
#define _elaspect_time_stepping_function_h

#include <elaspect/time_stepping/interface.h>
#include <elaspect/simulator_access.h>

#include <deal.II/base/parsed_function.h>

namespace elaspect
{
  namespace TimeStepping
  {
    template <int dim>
    class Function : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        void update () override;

        double
        get_next_time_step_size () const override;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        Functions::ParsedFunction<dim> time_stepping_function;
    };
  }
}

#endif
