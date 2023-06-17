#ifndef _elaspect_initial_topography_interface_h
#define _elaspect_initial_topography_interface_h

#include <elaspect/global.h>
#include <elaspect/plugins.h>

namespace elaspect
{
  namespace InitialTopography
  {
    using namespace dealii;

    template <int dim>
    class Interface
    {
      public:
        virtual ~Interface() = default;

        virtual void initialize ();

        virtual
        double value (const Point<dim-1> &p) const = 0;

        virtual
        double max_topography () const = 0;

        static
        void 
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);
    };

    template <int dim>
    void
    register_initial_topography_model (const std::string &name,
                                       const std::string &description,
                                       void (*declare_parameters_function) (ParameterHandler &),
                                       Interface<dim> *(*factory_function) ());

    template <int dim>
    Interface<dim> *
    create_initial_topography_model (ParameterHandler &prm);

    template <int dim>
    void
    declare_parameters (ParameterHandler &prm);

#define ELASPECT_REGISTER_INITIAL_TOPOGRAPHY_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ELASPECT_REGISTER_INITIAL_TOPOGRAPHY_MODEL_ ## classname \
  { \
    elaspect::internal::Plugins::RegisterHelper<elaspect::InitialTopography::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&elaspect::InitialTopography::register_initial_topography_model<2>, \
                                name, description); \
    elaspect::internal::Plugins::RegisterHelper<elaspect::InitialTopography::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&elaspect::InitialTopography::register_initial_topography_model<3>, \
                                name, description); \
  }
  }
}

#endif
