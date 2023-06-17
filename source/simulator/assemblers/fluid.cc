#include <elaspect/simulator/assemblers/hydro_system.h>

namespace elaspect
{
  namespace Assemblers
  {
    template <int dim>
    void
    HydroUnsaturatedFlow<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &/*scratch_base*/,
             internal::Assembly::CopyData::CopyDataBase<dim> &/*data_base*/) const
    {

    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace Assemblers
  {
#define INSTANTIATE(dim) \
    template class HydroUnsaturatedFlow<dim>;

    ELASPECT_INSTANTIATE(INSTANTIATE)
  }
}
