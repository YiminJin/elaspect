#ifndef _elaspect_boundary_velocity_interface_h
#define _elaspect_boundary_velocity_interface_h

#include <elaspect/plugins.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace BoundaryVelocity
  {
    template <int dim>
    class Interface
    {
      public:
        virtual ~Interface ();

        virtual void initialize ();

        virtual void update ();

        virtual
        Tensor<1,dim>
        boundary_velocity (const types::boundary_id boundary_indicator,
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

        Tensor<1,dim>
        boundary_velocity (const types::boundary_id boundary_indicator,
                           const Point<dim> &position) const;

        static
        void
        register_boundary_velocity (const std::string &name,
                                    const std::string &description,
                                    void (*declare_parameters_function) (ParameterHandler &),
                                    Interface<dim> *(*factory_function) ());

        const std::map<types::boundary_id, std::pair<std::string,std::vector<std::string> > > &
        get_active_boundary_velocity_names () const;

        const std::map<types::boundary_id,std::vector<std::unique_ptr<BoundaryVelocity::Interface<dim> > > > &
        get_active_boundary_velocity_conditions () const;

        const std::set<types::boundary_id> &
        get_zero_boundary_velocity_indicators () const;

        const std::set<types::boundary_id> &
        get_tangential_boundary_velocity_indicators () const;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm);

      private:
        std::map<types::boundary_id,std::vector<std::unique_ptr<BoundaryVelocity::Interface<dim> > > > boundary_velocity_objects;

        std::map<types::boundary_id, std::pair<std::string,std::vector<std::string> > > boundary_velocity_indicators;

        std::set<types::boundary_id> zero_velocity_boundary_indicators;

        std::set<types::boundary_id> tangential_velocity_boundary_indicators;
    };


#define ELASPECT_REGISTER_BOUNDARY_VELOCITY_MODEL(classname, name, description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ELASPECT_REGISTER_BOUNDARY_VELOCITY_MODEL_ ## classname \
  { \
    elaspect::internal::Plugins::RegisterHelper<elaspect::BoundaryVelocity::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&elaspect::BoundaryVelocity::Manager<2>::register_boundary_velocity, \
                                name, description); \
    elaspect::internal::Plugins::RegisterHelper<elaspect::BoundaryVelocity::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&elaspect::BoundaryVelocity::Manager<3>::register_boundary_velocity, \
                                name, description); \
  }
  }
}

#endif
