#ifndef _elaspect_mesh_refinement_composition_threshold_h
#define _elaspect_mesh_refinement_composition_threshold_h

#include <elaspect/mesh_refinement/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace MeshRefinement
  {
    template <int dim>
    class CompositionThreshold : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * After cells have been marked for coarsening/refinement, apply
         * additional criteria independent of the error estimate.
         */
        void
        tag_additional_cells () const override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * The thresholds that should be used for the individual
         * compositional fields.
         */
        std::vector<double> composition_thresholds;       
    };
  }
}

#endif
