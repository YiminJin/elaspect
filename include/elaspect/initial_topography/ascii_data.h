#ifndef _elaspect_initial_topography_ascii_data_h
#define _elaspect_initial_topography_ascii_data_h

#include <elaspect/initial_topography/interface.h>
#include <elaspect/structured_data.h>

namespace elaspect
{
  namespace InitialTopography
  {
    template <int dim>
    class AsciiData : public Utilities::AsciiDataBoundary<dim>, public Interface<dim>
    {
      public:
        AsciiData();

        void
        initialize() override;

        using Utilities::AsciiDataBoundary<dim>::initialize;

        double
        value(const Point<dim-1> &surface_point) const override;

        double max_topography() const override;

        static
        void
        declare_parameters(ParameterHandler &prm);

        void
        parse_parameters(ParameterHandler &prm) override;

      private:
        types::boundary_id surface_boundary_id;
    };
  }
}

#endif
