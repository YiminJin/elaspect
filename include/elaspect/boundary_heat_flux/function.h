#ifndef _elaspect_boundary_heat_flux_function_h
#define _elaspect_boundary_heat_flux_function_h

#include <elaspect/boundary_heat_flux/interface.h>
#include <elaspect/simulator_access.h>
#include <elaspect/utilities.h>

#include <deal.II/base/parsed_function.h>

namespace elaspect
{
  namespace BoundaryHeatFlux
  {
    using namespace dealii;

    template <int dim>
    class Function : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        Tensor<1,dim>
        heat_flux (const types::boundary_id boundary_id,
                   const Point<dim> &position,
                   const Tensor<1,dim> &normal_vector) const override;

        void update () override;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        Functions::ParsedFunction<dim> boundary_heat_flux_function;

        Utilities::Coordinates::CoordinateSystem coordinate_system;
    };
  }
}


#endif
