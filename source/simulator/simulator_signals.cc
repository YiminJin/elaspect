#include <elaspect/simulator_signals.h>

namespace elaspect
{
  namespace internals
  {
    namespace SimulatorSignals
    {
      std::list<std::function<void (elaspect::SimulatorSignals<2> &)>> connector_functions_2d;
      std::list<std::function<void (elaspect::SimulatorSignals<3> &)>> connector_functions_3d;

      static bool connector_functions_have_been_called = false;

      // add a user-provided connector to the list of connectors we keep
      void register_connector_function_2d(const std::function<void (elaspect::SimulatorSignals<2> &)> &connector)
      {
        Assert(!connector_functions_have_been_called,
               ExcMessage("Registration of signal connector happened after connection has already been called!"));
        connector_functions_2d.push_back (connector);       
      }


      void register_connector_function_3d(const std::function<void (elaspect::SimulatorSignals<3> &)> &connector)
      {
        Assert(!connector_functions_have_been_called,
               ExcMessage("Registration of signal connector happened after connection has already been called!"));
        connector_functions_3d.push_back (connector);
      }


      // call connectors to ensure that plugins get a change to register their slots
      template <>
      void call_connector_functions(elaspect::SimulatorSignals<2> &signals)
      {
        Assert(!connector_functions_have_been_called, 
               ExcInternalError());
        
        for (const auto &p : connector_functions_2d)
          p(signals);

        connector_functions_have_been_called = true;
      }


      template <>
      void call_connector_functions(elaspect::SimulatorSignals<3> &signals)
      {
        Assert(!connector_functions_have_been_called,
               ExcInternalError());

        for (const auto &p : connector_functions_3d)
          p(signals);

        connector_functions_have_been_called = true;
      }
    }
  }
}


// explicit instantiations
namespace elaspect
{
#define INSTANTIATE(dim) \
  template struct SimulatorSignals<dim>;


  ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
