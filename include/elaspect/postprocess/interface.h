#ifndef _elaspect_postprocess_interface_h
#define _elaspect_postprocess_interface_h

#include <elaspect/global.h>
#include <elaspect/plugins.h>
#include <elaspect/simulator_access.h>

#include <deal.II/base/table_handler.h>
#include <deal.II/base/parameter_handler.h>

namespace elaspect
{
  namespace Postprocess
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
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) = 0;

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);

        virtual
        std::list<std::string>
        required_other_postprocessors () const;
    };

    template <int dim>
    class Manager : public SimulatorAccess<dim>
    {
      public:
        std::list<std::pair<std::string,std::string> >
        execute (TableHandler &statistics);

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm);

        static
        void
        register_postprocessor (const std::string &name,
                                const std::string &description,
                                void (*declare_parameters_function) (ParameterHandler &),
                                Interface<dim> *(*factory_function) ());

      private:
        std::vector<std::unique_ptr<Interface<dim> > > postprocessors;
    };

#define ELASPECT_REGISTER_POSTPROCESSOR(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ELASPECT_REGISTER_POSTPROCESSOR_ ## classname \
  { \
    elaspect::internal::Plugins::RegisterHelper<elaspect::Postprocess::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&elaspect::Postprocess::Manager<2>::register_postprocessor, \
                                name, description); \
    elaspect::internal::Plugins::RegisterHelper<elaspect::Postprocess::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&elaspect::Postprocess::Manager<3>::register_postprocessor, \
                                name, description); \
  }
  }
}

#endif
