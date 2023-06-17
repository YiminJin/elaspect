#ifndef _elaspect_postprocess_stress_statistics_h
#define _elaspect_postprocess_stress_statistics_h

#include <elaspect/postprocess/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace Postprocess
  {
    template <int dim>
    class StressStatistics : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;
    };
  }
}

#endif
