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


#ifndef _elaspect_introspection_h
#define _elaspect_introspection_h

#include <elaspect/fe_variable_collection.h>
#include <elaspect/parameters.h>

namespace elaspect
{
  using namespace dealii;

  template <int dim>
  std::vector<VariableDeclaration<dim> >
  construct_default_variables (const Parameters<dim> &parameters);

  template <int dim>
  struct Introspection : public FEVariableCollection<dim>
  {
    public:
      Introspection (const std::vector<VariableDeclaration<dim> > &variables,
                     const Parameters<dim>                        &parameters);

      ~Introspection ();

      const unsigned int n_components;

      const unsigned int n_compositional_fields;

      const unsigned int n_quadrature_points;

      struct ComponentIndices
      {
        unsigned int displacement[dim];
        unsigned int temperature;
        std::vector<unsigned int> compositional_fields;
        std::vector<unsigned int> stress;
      };

      const ComponentIndices component_indices;

      const unsigned int n_blocks;

      struct BlockIndices
      {
        unsigned int displacement;
        unsigned int temperature;
        std::vector<unsigned int> compositional_fields;
        std::vector<unsigned int> stress;
      };

      const BlockIndices block_indices;

      struct Extractors
      {
        Extractors (const ComponentIndices &component_indices);

        const FEValuesExtractors::Vector displacement;
        const FEValuesExtractors::Scalar temperature;
        const std::vector<FEValuesExtractors::Scalar> compositional_fields;
        FEValuesExtractors::SymmetricTensor<2> stress;
      };

      const Extractors extractors;

      struct BaseElements
      {
        unsigned int displacement;
        unsigned int temperature;
        unsigned int compositional_fields;
        unsigned int stress;
      };

      const BaseElements base_elements;

      struct IndexSets
      {
        IndexSet system_relevant_set;

        std::vector<IndexSet> system_partitioning;

        std::vector<IndexSet> system_relevant_partitioning;
      };

      IndexSets index_sets;

      struct ComponentMasks
      {
        ComponentMasks (FEVariableCollection<dim> &fevs);

        ComponentMask displacement;
        ComponentMask temperature;
        std::vector<ComponentMask> compositional_fields;
        std::vector<ComponentMask> stress;
      };

      const ComponentMasks component_masks;

      struct QPDIndicators
      {
        QPDIndicators (const Parameters<dim> &parameters);

        unsigned int composition;
        unsigned int old_stress;
        unsigned int current_stress;
        unsigned int old_plastic_strain;
        unsigned int current_plastic_strain;
        unsigned int return_mapping_flag;
        unsigned int old_temperature;
      };

      const QPDIndicators qpd_indicators;

      const unsigned int n_qpd_components;

      std::vector<types::global_dof_index> system_dofs_per_block;

      unsigned int 
      compositional_index_for_name (const std::string &name) const;

      std::string
      name_for_compositional_index (const unsigned int index) const;

    private:
      std::vector<std::string> composition_names;
  };
}

#endif
