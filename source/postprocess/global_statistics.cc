#include <elaspect/postprocess/global_statistics.h>

namespace elaspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    GlobalStatistics<dim>::execute (TableHandler &/*statistics*/)
    {
      return std::make_pair (std::string(),std::string());
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace Postprocess
  {
    ELASPECT_REGISTER_POSTPROCESSOR(GlobalStatistics,
                                "global statistics",
                                "A postprocessor that outputs all the global statistics "
                                "information, e.g. the time of the simulation, the timestep "
                                "number, number of degrees of freedom and solver iterations "
                                "for each timestep. The postprocessor can output different "
                                "formats, the first printing one line in the statistics file "
                                "per nonlinear solver iteration (if a nonlinear solver scheme "
                                "is selected). The second prints one line per timestep, "
                                "summing the information about all nonlinear iterations in "
                                "this line. Note that this postprocessor is always active "
                                "independent on whether or not it is selected in the "
                                "parameter file.")
  }
}
