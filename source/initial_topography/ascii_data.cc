#include <elaspect/initial_topography/ascii_data.h>
#include <elaspect/geometry_model/box.h>
#include <elaspect/geometry_model/sphere.h>
#include <elaspect/geometry_model/spherical_shell.h>

namespace elaspect
{
  namespace InitialTopography
  {
    template <int dim>
    AsciiData<dim>::AsciiData()
      : surface_boundary_id(-1)
    {}


    template <int dim>
    void AsciiData<dim>::initialize()
    {
      surface_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");

      Utilities::AsciiDataBoundary<dim>::initialize({surface_boundary_id}, 1);
    }


    template <int dim>
    double
    AsciiData<dim>::value(const Point<dim-1> &surface_point) const
    {
      // Because the get_data_component function of AsciiDataBoundary
      // expects a dim-dimensional cartesian point, we have to
      // add a coordinate here, and, for spherical geometries,
      // change to cartesian coordinates.
      Point<dim> global_point;
      if (Plugins::plugin_type_matches<const GeometryModel::Box<dim>>(this->get_geometry_model()))
      {
        // No need to set the vertical coordinate correctly,
        // because it will be thrown away in get_data_component anyway
        for (unsigned int d = 0; d < dim-1; ++d)
          global_point[d] = surface_point[d];
      }
      else if (Plugins::plugin_type_matches<const GeometryModel::Sphere<dim>>(this->get_geometry_model()) ||
               Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>>(this->get_geometry_model()))
      {
        // No need to set the radial coordinate correctly,
        // because it will be thrown away in get_data_component anyway
        std::array<double, dim> point;
        for (unsigned int d = 0; d < dim-1; ++d)
          point[d+1] = surface_point[d];

        global_point = Utilities::Coordinates::spherical_to_cartesian_coordinates<dim>(point);
      }
      else
        AssertThrow(false, ExcNotImplemented());

      const double topo = Utilities::AsciiDataBoundary<dim>::get_data_component(surface_boundary_id,
                                                                                global_point,
                                                                                0);

      return topo;
    }


    template <int dim>
    double
    AsciiData<dim>::max_topography() const
    {
      return Utilities::AsciiDataBoundary<dim>::get_maximum_component_value(surface_boundary_id, 0);
    }


    template <int dim>
    void
    AsciiData<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Initial topography");
      {
        Utilities::AsciiDataBase<dim>::declare_parameters(prm, "", "");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    AsciiData<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Initial topography");
      {
        Utilities::AsciiDataBase<dim>::parse_parameters(prm);
      }
      prm.leave_subsection();
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace InitialTopography
  {
    ELASPECT_REGISTER_INITIAL_TOPOGRAPHY_MODEL(AsciiData,
                                           "ascii data",
                                           "")
  }
}
