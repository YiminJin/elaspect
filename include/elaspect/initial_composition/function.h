#ifndef _elaspect_initial_composition_function_h
#define _elaspect_initial_composition_function_h

#include <elaspect/initial_composition/interface.h>

#include <deal.II/base/parsed_function.h>

namespace elaspect
{
  namespace InitialComposition
  {
    using namespace dealii;

    template <int dim>
    class Function : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        double
        initial_composition (const Point<dim> &position,
                             const unsigned int n_comp) const override;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        std::unique_ptr<Functions::ParsedFunction<dim> > function;

        Utilities::Coordinates::CoordinateSystem coordinate_system;
    };
  }
}

#endif
