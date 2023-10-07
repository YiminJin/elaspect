/*
  Copyright (C) 2023 by Yimin Jin.

  This file is part of elASPECT.

  elASPECT is modified from the free software ASPECT; you can 
  redistribute it and/or modify it under the terms of the GNU 
  General Public License as published by the Free Software 
  Foundation; either version 2, or (at your option) any later 
  version.

  elASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with elASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef _elaspect_postprocess_visualization_material_properties_h
#define _elaspect_postprocess_visualization_material_properties_h

#include <elaspect/postprocess/visualization.h>
#include <elaspect/simulator_access.h>
#include <elaspect/material_model/io_interface.h>

#include <deal.II/numerics/data_postprocessor.h>

namespace elaspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      class MaterialProperties
        : public DataPostprocessor<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          MaterialProperties();

          std::vector<std::string>
          get_names() const override;

          std::vector<DataComponentInterpretation::DataComponentInterpretation>
          get_data_component_interpretation() const override;

          UpdateFlags
          get_needed_update_flags() const override;

          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double>> &computed_quantities) const override;

          static
          void
          declare_parameters(ParameterHandler &prm);

          void
          parse_parameters(ParameterHandler &prm) override;

        private:
          std::vector<std::string> property_names;

          MaterialModel::FieldDependences::Dependence field_dependences;
          MaterialModel::MaterialProperties::Property output_properties;
      };
    }
  }
}

#endif
