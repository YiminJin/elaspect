#ifndef _elaspect_postprocess_global_statistics_h
#define _elaspect_postprocess_global_statistics_h

#include <elaspect/postprocess/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace Postprocess
  {
    template <int dim>
    class GlobalStatistics : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;
    };
  }
}

#endif
