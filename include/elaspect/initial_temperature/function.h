#ifndef _elaspect_initial_temperature_function_h
#define _elaspect_initial_temperature_function_h

#include <elaspect/initial_temperature/interface.h>
#include <elaspect/simulator_access.h>

#include <deal.II/base/parsed_function.h>

namespace elaspect
{
  namespace InitialTemperature
  {
    using namespace dealii;

    template <int dim>
    class Function : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        Function ();

        double initial_temperature (const Point<dim> &position) const override;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        Functions::ParsedFunction<dim> function;

        Utilities::Coordinates::CoordinateSystem coordinate_system;
    };
  }
}

#endif
