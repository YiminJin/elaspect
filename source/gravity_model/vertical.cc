#include <elaspect/gravity_model/vertical.h>
#include <elaspect/geometry_model/interface.h>

#include <deal.II/base/tensor.h>

namespace elaspect
{
  namespace GravityModel
  {
    template <int dim>
    Tensor<1,dim>
    Vertical<dim>::gravity_vector (const Point<dim> &) const
    {
      Tensor<1,dim> g;
      g[dim-1] = -gravity_magnitude;
      return g;
    }
    template <int dim>
    void
    Vertical<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Vertical");
        {
          prm.declare_entry ("Magnitude", "1.",
                             Patterns::Double (),
                             "Value of the gravity vector in $m/s^2$ directed "
                             "along negative y (2D) or z (3D) axis (if the magnitude "
                             "is positive.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    Vertical<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Vertical");
        {
          gravity_magnitude = prm.get_double ("Magnitude");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

      AssertThrow (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::cartesian,
                   ExcMessage ("Gravity model 'vertical' should not be used with geometry models that "
                               "do not have a cartesian natural coordinate system."));
    }
  }
}

// explicit instantiations
namespace elaspect
{
  namespace GravityModel
  {
    ELASPECT_REGISTER_GRAVITY_MODEL(Vertical,
                                "vertical",
                                "A gravity model in which the gravity direction is vertical (pointing "
                                "downward for positive values) and at a constant magnitude by default "
                                "equal to one.")
  }
}
