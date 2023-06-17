#ifndef _elaspect_termination_criteria_end_time_h
#define _elaspect_termination_criteria_end_time_h

#include <elaspect/termination_criteria/interface.h>

namespace elaspect
{
  namespace TerminationCriteria
  {
    template <int dim>
    class EndTime : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        double 
        check_for_last_time_step (const double time_step) const override;
        
        bool
        execute () override;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        double end_time;
    };
  }
}

#endif
