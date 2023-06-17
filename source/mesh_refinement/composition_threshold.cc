#include <elaspect/mesh_refinement/composition_threshold.h>
#include <elaspect/quadrature_point_data.h>

namespace elaspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    CompositionThreshold<dim>::tag_additional_cells() const
    {
      AssertThrow (this->n_compositional_fields() >= 1,
                   ExcMessage ("Refinement criterion 'composition threshold' can not "
                               "be used when no compositional fields are active!"));

      std::vector<double> composition_values(this->get_qpd_handler().n_quadrature_points());

      const Introspection<dim> &introspection = this->introspection();

      for (const auto &cell : this->get_qpd_handler().active_cell_iterators())
        if (cell->is_locally_owned())
        {
          bool refine = false;

          for (unsigned int c = 0; c < introspection.n_compositional_fields; ++c)
          {
            cell->get(introspection.qpd_indicators.composition + c, composition_values);
            for (unsigned int q = 0; q < composition_values.size(); ++q)
              if (composition_values[q] > composition_thresholds[c])
              {
                refine = true;
                break;
              }
          }

          if (refine)
          {
            cell->clear_coarsen_flag();
            cell->set_refine_flag();
          }
        }
    }


    template <int dim>
    void
    CompositionThreshold<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Composition threshold");
        {
          prm.declare_entry("Compositional field thresholds", "",
                            Patterns::List (Patterns::Double()),
                            "A list of thresholds that every individual compositional "
                            "field will be evaluated against.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    CompositionThreshold<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Composition threshold");
        {
          composition_thresholds = Utilities::string_to_double(
                Utilities::split_string_list(prm.get("Compositional field thresholds")));

          AssertThrow (composition_thresholds.size() == this->n_compositional_fields(),
                       ExcMessage ("The number of thresholds given here must be "
                                   "equal to the number of chosen refinement criteria."));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace MeshRefinement
  {
    ELASPECT_REGISTER_MESH_REFINEMENT_CRITERION(CompositionThreshold,
                                            "composition threshold",
                                            "A mesh refinement criterion that computes refinement "
                                            "indicators from the compositional fields. If any field "
                                            "exceeds the threshold given in the input file, the cell "
                                            "is marked for refinement.")
  }
}
