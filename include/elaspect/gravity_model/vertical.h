#ifndef _elaspect_gravity_model_vertical_h
#define _elaspect_gravity_model_vertical_h

#include <elaspect/gravity_model/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace GravityModel
  {
    using namespace dealii;

    template <int dim>
    class Vertical : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the gravity vector as a function of position.
         */
        Tensor<1,dim> gravity_vector (const Point<dim> &position) const override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * Magnitude of the gravity vector.
         */
        double gravity_magnitude;

    };
  }
}

#endif
