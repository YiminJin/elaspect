#ifndef _elaspect_termination_criteria_end_step_h
#define _elaspect_termination_criteria_end_step_h

#include <elaspect/termination_criteria/interface.h>

namespace elaspect
{
  namespace TerminationCriteria
  {
    template <int dim>
    class EndStep : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        bool
        execute () override;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        unsigned int end_step;
    };
  }
}

#endif
