#ifndef _elaspect_gravity_model_radial_h
#define _elaspect_gravity_model_radial_h

#include <elaspect/gravity_model/interface.h>
#include <elaspect/simulator_access.h>
#include <elaspect/utilities.h>

#include <deal.II/base/parsed_function.h>

namespace elaspect
{
  namespace GravityModel
  {
    template <int dim>
    class RadialConstant : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        Tensor<1,dim> gravity_vector (const Point<dim> &position) const override;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        double magnitude;
    };


    template <int dim>
    class RadialLinear : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        Tensor<1,dim> gravity_vector (const Point<dim> &position) const override;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        double magnitude_at_surface;

        double magnitude_at_bottom;
    };

    template <int dim>
    class RadialFunction : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        Tensor<1,dim> gravity_vector (const Point<dim> &position) const override;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        Functions::ParsedFunction<1> magnitude_function;
    };
  }
}

#endif
