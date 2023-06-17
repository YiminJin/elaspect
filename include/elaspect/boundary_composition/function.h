#ifndef _elaspect_boundary_composition_function_h
#define _elaspect_boundary_composition_function_h

#include <elaspect/boundary_composition/interface.h>
#include <elaspect/simulator_access.h>
#include <elaspect/utilities.h>

#include <deal.II/base/parsed_function.h>

namespace elaspect
{
  namespace BoundaryComposition
  {
    using namespace dealii;

    template <int dim>
    class Function : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        double boundary_composition (const types::boundary_id boundary_indicator,
                                     const Point<dim> &position,
                                     const unsigned int compositional_field) const override;

        void update () override;

        static
        void declare_parameters (ParameterHandler &prm);

        void parse_parameters (ParameterHandler &prm) override;

      private:
        std::unique_ptr<Functions::ParsedFunction<dim>> function;

        Utilities::Coordinates::CoordinateSystem coordinate_system;
    };
  }
}

#endif
