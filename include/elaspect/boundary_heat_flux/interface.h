#ifndef _elaspect_boundary_heat_flux_interface_h
#define _elaspect_boundary_heat_flux_interface_h

#include <elaspect/global.h>
#include <elaspect/plugins.h>

namespace elaspect
{
  namespace BoundaryHeatFlux
  {
    using namespace dealii;

    template <int dim>
    class Interface
    {
      public:
        virtual ~Interface ();

        virtual void initialize ();

        virtual void update ();

        virtual
        Tensor<1,dim>
        heat_flux (const types::boundary_id boundary_id,
                   const Point<dim> &position,
                   const Tensor<1,dim> &normal_vector) const = 0;

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);
    };

    template <int dim>
    void
    register_boundary_heat_flux (const std::string &name,
                                 const std::string &description,
                                 void (*declare_parameters_function) (ParameterHandler &),
                                 Interface<dim> *(*factory_function) ());

    template <int dim>
    Interface<dim> *
    create_boundary_heat_flux (ParameterHandler &prm);

    template <int dim>
    void
    declare_parameters (ParameterHandler &prm);

#define ELASPECT_REGISTER_BOUNDARY_HEAT_FLUX_MODEL(classname, name, description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ELASPECT_REGISTER_BOUNDARY_HEAT_FLUX_MODEL_ ## classname \
  { \
    elaspect::internal::Plugins::RegisterHelper<elaspect::BoundaryHeatFlux::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&elaspect::BoundaryHeatFlux::register_boundary_heat_flux<2>, \
                                name, description); \
    elaspect::internal::Plugins::RegisterHelper<elaspect::BoundaryHeatFlux::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&elaspect::BoundaryHeatFlux::register_boundary_heat_flux<3>, \
                                name, description); \
  }
  }
}

#endif
