#include <elaspect/postprocess/visualization/material_properties.h>
#include <elaspect/material_model/handler.h>
#include <elaspect/utilities.h>

namespace elaspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      MaterialProperties<dim>::MaterialProperties()
        : DataPostprocessor<dim>()
      {}


      template <int dim>
      std::vector<std::string>
      MaterialProperties<dim>::get_names() const
      {
        std::vector<std::string> solution_names;

        for (const auto &property_name : property_names)
        {
          if (property_name == "reaction terms")
          {
            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              solution_names.push_back(this->introspection().name_for_compositional_index(c) + "_change");
          }
          else
          {
            solution_names.push_back(property_name);
            std::replace(solution_names.back().begin(), solution_names.back().end(), ' ', '_');
          }
        }

        return solution_names;
      }


      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      MaterialProperties<dim>::get_data_component_interpretation () const
      {
        std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation;
        for (const auto &property_name : property_names)
        {
          if (property_name == "reaction terms")
          {
            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              interpretation.push_back (DataComponentInterpretation::component_is_scalar);
          }
          else
            interpretation.push_back (DataComponentInterpretation::component_is_scalar);
        }

        return interpretation;
      }


      template <int dim>
      UpdateFlags
      MaterialProperties<dim>::get_needed_update_flags () const
      {
        return update_gradients | update_values  | update_quadrature_points;
      }



      template <int dim>
      void
      MaterialProperties<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert(computed_quantities.size() == n_quadrature_points,
               ExcInternalError());
        Assert(input_data.solution_values[0].size() == this->introspection().n_components,
               ExcInternalError());

        MaterialModel::MaterialModelInputs<dim> in(input_data,
                                                   this->introspection(),
                                                   field_dependences,
                                                   output_properties);
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points,
                                                     this->n_compositional_fields());

        this->get_material_handler().evaluate(in, out);

        // We want to output material properties as they are used in the
        // program during assembly. To do so, some of the material averaging
        // modes require a quadrature object -- but we do not have this,
        // all we have is the mapped quadrature points. As a consequence,
        // only do the averaging for those modes that do not require a
        // functional quadrature object. This means that we only
        // do the wrong thing for the Q1 averaging where the difference
        // between input and output of the averaging operation is generally
        // small and probably not visible anyway.
        //
        // The average() function checks whether any quadrature object
        // it uses has the correct size. So passing an invalid object
        // in the following code carries little risk if the list of
        // averaging modes that require a quadrature object expands:
        // Every time we generate graphical output for a model that
        // uses this kind of averaging, we will trigger an exception.
        if (this->get_parameters().material_averaging != 
            MaterialModel::MaterialAveraging::AveragingOperation::project_to_Q1)
          MaterialModel::MaterialAveraging::average(this->get_parameters().material_averaging,
                                                    input_data.template get_cell<dim>(),
                                                    Quadrature<dim>(),
                                                    this->get_mapping(),
                                                    false,
                                                    out);

        for (unsigned int q = 0; q < n_quadrature_points; ++q)
        {
          unsigned int output_index = 0;
          for (unsigned int i = 0; i < property_names.size(); ++i, ++output_index)
          {
            if (property_names[i] == "viscosity")
              computed_quantities[q][output_index] = out.viscosities[q];

            else if (property_names[i] == "density")
              computed_quantities[q][output_index] = out.densities[q];

            else if (property_names[i] == "specific heat")
              computed_quantities[q][output_index] = out.specific_heat[q];

            else if (property_names[i] == "thermal expansivity")
              computed_quantities[q][output_index] = out.thermal_expansivities[q];

            else if (property_names[i] == "thermal conductivity")
              computed_quantities[q][output_index] = out.thermal_conductivities[q];

            else if (property_names[i] == "thermal diffusivity")
              computed_quantities[q][output_index] = out.thermal_conductivities[q] /
                                                     (out.densities[q] * out.specific_heat[q]);
            
            else if (property_names[i] == "entropy derivative pressure")
              computed_quantities[q][output_index] = out.entropy_derivative_pressure[q];

            else if (property_names[i] == "entropy derivative temperature")
              computed_quantities[q][output_index] = out.entropy_derivative_temperature[q];

            else if (property_names[i] == "reaction terms")
            {
              for (unsigned int k = 0; k < this->n_compositional_fields(); ++k, ++output_index)
                computed_quantities[q][output_index] = out.reaction_terms[q][k];

              --output_index;
            }
            else
              AssertThrow(false,
                          ExcMessage("Material property <" + property_names[i] +
                                     "< not implemented for visualization."));
          }
        }
      }


      template <int dim>
      void
      MaterialProperties<dim>::declare_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Material properties");
            {
              const std::string pattern_of_names
                = "viscosity|density|specific heat|thermal expansivity|"
                  "thermal conductivity|thermal diffusivity|"
                  "entropy derivative pressure|entropy derivative temperature|"
                  "reaction terms";

              prm.declare_entry("List of material properties",
                                "density,thermal expansivity,specific heat,viscosity",
                                Patterns::MultipleSelection(pattern_of_names),
                                "A comma separated list of material properties that should be "
                                "written whenever writing graphical output. "
                                "The following material properties are available:\n\n"
                                +
                                pattern_of_names);             
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      MaterialProperties<dim>::parse_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Material properties");
            {
              // Get property names and compare against variable names
              property_names = Utilities::split_string_list(prm.get("List of material properties"));
              AssertThrow(Utilities::has_unique_entries(property_names),
                          ExcMessage("The list of strings for the parameter "
                                     "'Postprocess/Visualization/Material properties/List of material properties' "
                                     "contains entries more than once. This is not allowed. "
                                     "Please check your parameter file."));

              // Get the field dependences
              output_properties = MaterialModel::MaterialProperties::none;
              for (unsigned int i = 0; i < property_names.size(); ++i)
              {
                if (property_names[i] == "viscosity")
                  output_properties |= MaterialModel::MaterialProperties::viscosity;

                else if (property_names[i] == "density")
                  output_properties |= MaterialModel::MaterialProperties::density;

                else if (property_names[i] == "specific heat")
                  output_properties |= MaterialModel::MaterialProperties::specific_heat;

                else if (property_names[i] == "thermal expansivity")
                  output_properties |= MaterialModel::MaterialProperties::thermal_expansivity;

                else if (property_names[i] == "thermal conductivity" ||
                         property_names[i] == "thermal diffusivity")
                  output_properties |= MaterialModel::MaterialProperties::thermal_conductivity;

                else if (property_names[i] == "entropy derivative pressure")
                  output_properties |= MaterialModel::MaterialProperties::entropy_derivative_pressure;

                else if (property_names[i] == "entropy derivative temperature")
                  output_properties |= MaterialModel::MaterialProperties::entropy_derivative_temperature;

                else if (property_names[i] == "reaction terms")
                  output_properties |= MaterialModel::MaterialProperties::reaction_terms;

                else 
                  AssertThrow(false,
                              ExcMessage("Material property <" + property_names[i] +
                                         "< not implemented for visualization."));
              }

              field_dependences = this->get_material_handler().get_field_dependences_for_evaluation(output_properties);
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ELASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(MaterialProperties,
                                                "material properties",
                                                "")
    }
  }
}
