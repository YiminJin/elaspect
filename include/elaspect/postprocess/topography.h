#ifndef _elaspect_postprocess_topography_h
#define _elaspect_postprocess_topography_h

#include <elaspect/postprocess/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace Postprocess
  {
    template <int dim>
    class Topography : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        std::pair<std::string,std::string> execute (TableHandler &statistics) override;

        static
        void
        declare_parameters(ParameterHandler &prm);

        void
        parse_parameters(ParameterHandler &prm) override;

      private:
        bool write_to_file;

        double output_interval;

        double last_output_time;
    };
  }
}

#endif
