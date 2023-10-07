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


#include <elaspect/postprocess/visualization/partition.h>

namespace elaspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      Partition<dim>::Partition()
        :
        DataPostprocessorScalar<dim>("partition",
                                     update_default)
      {}


      template <int dim>
      void
      Partition<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &/*input_data*/,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        Assert(computed_quantities[0].size() == 1, ExcInternalError());

        for (auto &quantity : computed_quantities)
        {
          // simply get the partition number from the triangulation
          quantity(0) = this->get_triangulation().locally_owned_subdomain();
        }
      }
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ELASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(Partition,
                                                "partition",
                                                 "A visualization output object that generates output "
                                                 "for the parallel partition that every cell of the "
                                                 "mesh is associated with.")
    }
  }
}
