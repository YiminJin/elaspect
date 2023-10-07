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


#ifndef _elaspect_material_model_handler_h
#define _elaspect_material_model_handler_h

#include <elaspect/parameters.h>
#include <elaspect/material_model/io_interface.h>
#include <elaspect/material_model/utilities.h>
#include <elaspect/material_model/basic_properties/interface.h>
#include <elaspect/material_model/viscosity/interface.h>
#include <elaspect/material_model/plasticity/interface.h>

#include <deal.II/base/symmetric_tensor.h>

namespace elaspect
{
  template <int dim> class SimulatorAccess;

  using namespace dealii;

  template <int dim>
  class MaterialHandler : public SimulatorAccess<dim>
  {
    public:
      void initialize ();

      void update ();

      void 
      evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
               MaterialModel::MaterialModelOutputs<dim> &out) const;

      void
      apply_return_mapping(const MaterialModel::MaterialModelInputs<dim> &in,
                           std::vector<SymmetricTensor<2,dim>> &          stresses,
                           std::vector<double> &                          plastic_strains,
                           std::vector<double> &                          flags) const;

      MaterialModel::FieldDependences::Dependence
      get_field_dependences_for_evaluation(const MaterialModel::MaterialProperties::Property properties) const;

      MaterialModel::FieldDependences::Dependence
      get_field_dependences_for_return_mapping() const;

      MaterialModel::MaterialProperties::Property
      get_material_properties_for_return_mapping() const;

      static
      void
      declare_parameters(ParameterHandler &prm);

      void
      parse_parameters(ParameterHandler &prm);

    private:
      MaterialUtilities::CompositionalAveragingOperation rheological_parameters_averaging;

      std::unique_ptr<MaterialModel::BasicProperties::Interface<dim>> basic_property_model;
      std::unique_ptr<MaterialModel::Plasticity::Interface<dim>>      plasticity_model;
      std::unique_ptr<MaterialModel::Viscosity::Interface<dim>>       viscosity_model;
  };
}

#endif
