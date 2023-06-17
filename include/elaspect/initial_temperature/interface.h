#ifndef _elaspect_initial_temperature_interface_h
#define _elaspect_initial_temperature_interface_h

#include <elaspect/plugins.h>
#include <elaspect/utilities.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace InitialTemperature
  {
    using namespace dealii;

    template <int dim>
    class Interface
    {
      public:
        virtual ~Interface ();

        virtual void initialize ();

        virtual
        double
        initial_temperature (const Point<dim> &position) const = 0;

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

        double
        initial_temperature (const Point<dim> &position) const;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm);

        static
        void
        register_initial_temperature (const std::string &name,
                                      const std::string &description,
                                      void (*declare_parameters_function) (ParameterHandler &),
                                      Interface<dim> *(*factory_function) ());

        const std::vector<std::string> &
        get_active_initial_temperature_names () const;

        const std::list<std::unique_ptr<Interface<dim> > > &
        get_active_initial_temperature_conditions () const;

      private:
        std::list<std::unique_ptr<Interface<dim> > > initial_temperature_objects;

        std::vector<std::string> model_names;

        std::vector<elaspect::Utilities::Operator> model_operators;
    };

    template <int dim>
    std::string
    get_valid_model_names_pattern ();


    /**
     * Given a class name, a name, and a description for the parameter file
     * for a initial conditions model, register it with the functions that can
     * declare their parameters and create these objects.
     *
     * @ingroup InitialTemperatures
     */
#define ELASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ELASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL_ ## classname \
  { \
    elaspect::internal::Plugins::RegisterHelper<elaspect::InitialTemperature::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&elaspect::InitialTemperature::Manager<2>::register_initial_temperature, \
                                name, description); \
    elaspect::internal::Plugins::RegisterHelper<elaspect::InitialTemperature::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&elaspect::InitialTemperature::Manager<3>::register_initial_temperature, \
                                name, description); \
  }
  }
}

#endif
