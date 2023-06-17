#ifndef _elaspect_termination_criteria_interface_h
#define _elaspect_termination_criteria_interface_h

#include <elaspect/plugins.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace TerminationCriteria
  {
    template <int dim>
    class Interface
    {
      public:
        virtual ~Interface ();

        virtual void initialize ();

        virtual bool execute () = 0;

        virtual double check_for_last_time_step (const double time_step) const;

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
        bool execute () const;

        double check_for_last_time_step (const double time_step) const;

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);

        static
        void
        register_termination_criterion (const std::string &name,
                                        const std::string &description,
                                        void (*declare_parameters_function) (ParameterHandler &),
                                        Interface<dim> *(*factory_function) ());

      private:
        std::list<std::unique_ptr<Interface<dim> > > termination_objects;

        std::list<std::string>                       termination_obj_names;
    };


    /**
     * Given a class name, a name, and a description for the parameter file
     * for a termination criterion object, register it with the
     * elaspect::TerminationCriteria::Manager class.
     *
     * @ingroup TerminationCriteria
     */
#define ELASPECT_REGISTER_TERMINATION_CRITERION(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ELASPECT_REGISTER_TERMINATION_CRITERION_ ## classname \
  { \
    elaspect::internal::Plugins::RegisterHelper<elaspect::TerminationCriteria::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&elaspect::TerminationCriteria::Manager<2>::register_termination_criterion, \
                                name, description); \
    elaspect::internal::Plugins::RegisterHelper<elaspect::TerminationCriteria::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&elaspect::TerminationCriteria::Manager<3>::register_termination_criterion, \
                                name, description); \
  }
  }
}

#endif
