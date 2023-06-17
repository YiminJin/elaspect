#ifndef _elaspect_boundary_composition_interface_h
#define _elaspect_boundary_composition_interface_h

#include <elaspect/plugins.h>
#include <elaspect/utilities.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace BoundaryComposition
  {
    using namespace dealii;

    template <int dim>
    class Interface
    {
      public:
        virtual ~Interface() = default;

        virtual void initialize();

        virtual void update();

        virtual
        double 
        boundary_composition(const types::boundary_id boundary_indicator,
                             const Point<dim> &position,
                             const unsigned int compositional_field) const = 0;

        static
        void
        declare_parameters(ParameterHandler &prm);

        virtual
        void
        parse_parameters(ParameterHandler &prm);
    };

    template <int dim>
    class Manager : public SimulatorAccess<dim>
    {
      public:
        ~Manager() override;

        virtual
        void
        update();

        static
        void
        declare_parameters(ParameterHandler &prm);

        void
        parse_parameters(ParameterHandler &prm);

        double
        boundary_composition(const types::boundary_id boundary_indicator,
                             const Point<dim> &position,
                             const unsigned int compositional_field) const;

        static
        void
        register_boundary_composition(const std::string &name,
                                      const std::string &description,
                                      void (*declare_parameters_function)(ParameterHandler &),
                                      Interface<dim> *(*factory_function)());

        const std::set<types::boundary_id> &
        get_fixed_composition_boundary_indicators() const;

      private:
        std::vector<std::unique_ptr<Interface<dim>>> boundary_composition_objects;

        std::vector<std::string> model_names;

        std::vector<Utilities::Operator> model_operators;

        std::set<types::boundary_id> fixed_composition_boundary_indicators;
    };

    template <int dim>
    std::string
    get_valid_model_names_pattern();


    /**
     * Given a class name, a name, and a description for the parameter file
     * for a boundary composition model, register it with the functions that
     * can declare their parameters and create these objects.
     *
     * @ingroup BoundaryCompositions
     */
#define ELASPECT_REGISTER_BOUNDARY_COMPOSITION_MODEL(classname, name, description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ELASPECT_REGISTER_BOUNDARY_COMPOSITION_MODEL_ ## classname \
  { \
    elaspect::internal::Plugins::RegisterHelper<elaspect::BoundaryComposition::Interface<2>,classname<2>> \
    dummy_ ## classname ## _2d (&elaspect::BoundaryComposition::Manager<2>::register_boundary_composition, \
                                name, description); \
    elaspect::internal::Plugins::RegisterHelper<elaspect::BoundaryComposition::Interface<3>,classname<3>> \
    dummy_ ## classname ## _3d (&elaspect::BoundaryComposition::Manager<3>::register_boundary_composition, \
                                name, description); \
  }
  }
}

#endif
