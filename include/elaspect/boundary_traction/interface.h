#ifndef _elaspect_boundary_traction_interface_h
#define _elaspect_boundary_traction_interface_h

#include <elaspect/plugins.h>
#include <elaspect/geometry_model/interface.h>

namespace elaspect
{
  namespace BoundaryTraction
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
        boundary_traction (const types::boundary_id boundary_indicator,
                           const Point<dim> &position,
                           const Tensor<1,dim> &normal_vector) const = 0;

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);

      protected:
        const GeometryModel::Interface<dim> *geometry_model;
    };

    template <int dim>
    void
    register_boundary_traction (const std::string &name,
                                const std::string &description,
                                void (*declare_parameters_function) (ParameterHandler &),
                                Interface<dim> *(*factory_function) ());

    template <int dim>
    Interface<dim> *
    create_boundary_traction (const std::string &name);

    template <int dim>
    std::string
    get_names ();

    template <int dim>
    void
    declare_parameters (ParameterHandler &prm);


    /**
     * Given a class name, a name, and a description for the parameter file
     * for a traction boundary conditions model, register it with the
     * functions that can declare their parameters and create these objects.
     *
     * @ingroup BoundaryTractions
     */
#define ELASPECT_REGISTER_BOUNDARY_TRACTION_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ELASPECT_REGISTER_BOUNDARY_TRACTION_MODEL_ ## classname \
  { \
    elaspect::internal::Plugins::RegisterHelper<elaspect::BoundaryTraction::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&elaspect::BoundaryTraction::register_boundary_traction<2>, \
                                name, description); \
    elaspect::internal::Plugins::RegisterHelper<elaspect::BoundaryTraction::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&elaspect::BoundaryTraction::register_boundary_traction<3>, \
                                name, description); \
  }
  }
}

#endif
