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
