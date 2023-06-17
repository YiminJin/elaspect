#ifndef _elaspect_mesh_refinement_interface_h
#define _elaspect_mesh_refinement_interface_h

#include <elaspect/plugins.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace MeshRefinement
  {
    using namespace dealii;

    template <int dim>
    class Interface
    {
      public:
        virtual
        ~Interface ();

        virtual void initialize ();

        virtual void update ();

        virtual
        void
        execute (Vector<float> &error_indicators) const;

        virtual
        void
        tag_additional_cells () const;

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
        void 
        execute (Vector<float> &error_indicators) const;

        virtual
        void
        tag_additional_cells () const;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm);

        static
        void
        register_mesh_refinement_criterion (const std::string &name,
                                            const std::string &description,
                                            void (*declare_parameters_function) (ParameterHandler &),
                                            Interface<dim> *(*factory_function) ());

      private:
        enum MergeOperation
        { plus, max };

        MergeOperation merge_operation;

        bool normalize_criteria;

        std::vector<double> scaling_factors;

        std::list<std::unique_ptr<Interface<dim> > > mesh_refinement_objects;
    };


#define ELASPECT_REGISTER_MESH_REFINEMENT_CRITERION(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ELASPECT_REGISTER_MESH_REFINEMENT_CRITERION_ ## classname \
  { \
    elaspect::internal::Plugins::RegisterHelper<elaspect::MeshRefinement::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&elaspect::MeshRefinement::Manager<2>::register_mesh_refinement_criterion, \
                                name, description); \
    elaspect::internal::Plugins::RegisterHelper<elaspect::MeshRefinement::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&elaspect::MeshRefinement::Manager<3>::register_mesh_refinement_criterion, \
                                name, description); \
  }
  }
}

#endif
