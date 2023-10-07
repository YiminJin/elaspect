/*
  Copyright (C) 2023 by Yimin Jin.

  This file is part of elASPECT.

  elASPECT is modified from the free software ASPECT; you can 
  redistribute it and/or modify it under the terms of the GNU 
  General Public License as published by the Free Software 
  Foundation; either version 2, or (at your option) any later 
  version.

  elASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with elASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


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
