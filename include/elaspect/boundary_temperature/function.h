#ifndef _elaspect_boundary_temperature_function_h
#define _elaspect_boundary_temperature_function_h

#include <elaspect/boundary_temperature/interface.h>
#include <elaspect/simulator_access.h>
#include <elaspect/utilities.h>

#include <deal.II/base/parsed_function.h>

namespace elaspect
{
  namespace BoundaryTemperature
  {
    using namespace dealii;

    template <int dim>
    class Function : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        Function ();

        void update () override;

        double
        boundary_temperature (const types::boundary_id boundary_indicator,
                              const Point<dim> &position) const override;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        Functions::ParsedFunction<dim> boundary_temperature_function;

        Utilities::Coordinates::CoordinateSystem coordinate_system;
    };
  }
}

#endif
