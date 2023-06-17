#ifndef _elaspect_boundary_temperature_interface_h
#define _elaspect_boundary_temperature_interface_h

#include <elaspect/plugins.h>
#include <elaspect/utilities.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace BoundaryTemperature
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
        double boundary_temperature (const types::boundary_id boundary_indicators,
                                     const Point<dim> &position) const = 0;

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);
    };


    template <int dim>
    class Manager : public SimulatorAccess<dim>
    {
      public:
        ~Manager () override;

        virtual void update ();

        virtual
        double boundary_temperature (const types::boundary_id boundary_indicators,
                                     const Point<dim> &position) const;

        static
        void
        register_boundary_temperature (const std::string &name,
                                       const std::string &description,
                                       void (*declare_parameters_function) (ParameterHandler &),
                                       Interface<dim> *(*factory_function) ());

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);

        const std::vector<std::string> &
        get_active_boundary_temperature_names () const;

        const std::vector<std::unique_ptr<Interface<dim> > > &
        get_active_boundary_temperature_conditions () const;

        const std::set<types::boundary_id> &
        get_fixed_temperature_boundary_indicators() const;

        bool
        allows_fixed_temperature_on_outflow_boundaries() const;

      private:
        std::vector<std::unique_ptr<Interface<dim> > > boundary_temperature_objects;

        std::vector<std::string> model_names;

        std::vector<elaspect::Utilities::Operator> model_operators;

        std::set<types::boundary_id> fixed_temperature_boundary_indicators;

        bool allow_fixed_temperature_on_outflow_boundaries;
    };

#define ELASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL(classname, name, description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ELASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL_ ## classname \
  { \
    elaspect::internal::Plugins::RegisterHelper<elaspect::BoundaryTemperature::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&elaspect::BoundaryTemperature::Manager<2>::register_boundary_temperature, \
                                name, description); \
    elaspect::internal::Plugins::RegisterHelper<elaspect::BoundaryTemperature::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&elaspect::BoundaryTemperature::Manager<3>::register_boundary_temperature, \
                                name, description); \
  }
  }
}

#endif
