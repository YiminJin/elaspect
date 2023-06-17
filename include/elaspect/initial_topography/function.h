#ifndef _elaspect_initial_topography_function_h
#define _elaspect_initial_topography_function_h

#include <elaspect/initial_topography/interface.h>
#include <elaspect/simulator_access.h>
#include <elaspect/utilities.h>

#include <deal.II/base/parsed_function.h>

namespace elaspect
{
  namespace InitialTopography
  {
    using namespace dealii;

    template <int dim>
    class Function : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        Function();

        double
        value(const Point<dim-1> &p) const override;

        double max_topography() const override;

        static
        void
        declare_parameters(ParameterHandler &prm);

        void 
        parse_parameters(ParameterHandler &prm) override;

      private:
        double max_topo;

        Functions::ParsedFunction<dim> initial_topography_function;

        Utilities::Coordinates::CoordinateSystem coordinate_system;
    };
  }
}

#endif
