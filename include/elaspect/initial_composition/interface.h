#ifndef _elaspect_initial_composition_interface_h
#define _elaspect_initial_composition_interface_h

#include <elaspect/plugins.h>
#include <elaspect/utilities.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace InitialComposition
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
        initial_composition (const Point<dim> &position,
                             const unsigned int n_comp) const = 0;

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
        initial_composition (const Point<dim> &position,
                             const unsigned int n_comp) const;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm);

        static
        void
        register_initial_composition (const std::string &name,
                                      const std::string &description,
                                      void (*declare_parameters_function) (ParameterHandler &),
                                      Interface<dim> *(*factory_function) ());

        const std::vector<std::string> &
        get_active_initial_composition_names () const;

        const std::list<std::unique_ptr<Interface<dim> > > &
        get_active_initial_composition_conditions () const;

      private:
        std::list<std::unique_ptr<Interface<dim> > > initial_composition_objects;

        std::vector<std::string> model_names;

        std::vector<elaspect::Utilities::Operator> model_operators;
    };


    /**
     * Given a class name, a name, and a description for the parameter file
     * for a initial composition model, register it with the functions that can
     * declare their parameters and create these objects.
     *
     * @ingroup InitialCompositions
     */
#define ELASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ELASPECT_REGISTER_INITIAL_COMPOSITION_MODEL_ ## classname \
  { \
    elaspect::internal::Plugins::RegisterHelper<elaspect::InitialComposition::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&elaspect::InitialComposition::Manager<2>::register_initial_composition, \
                                name, description); \
    elaspect::internal::Plugins::RegisterHelper<elaspect::InitialComposition::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&elaspect::InitialComposition::Manager<3>::register_initial_composition, \
                                name, description); \
  }
  }
}

#endif
