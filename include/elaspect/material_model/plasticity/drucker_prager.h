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


#ifndef _elaspect_material_model_plasticity_drucker_prager_h
#define _elaspect_material_model_plasticity_drucker_prager_h

#include <elaspect/material_model/plasticity/interface.h>
#include <elaspect/simulator_access.h>

#include <deal.II/base/parsed_function.h>

namespace elaspect
{
  namespace MaterialModel
  {
    namespace Plasticity
    {
      using namespace dealii;

      template <int dim>
      class DruckerPrager : public Interface<dim>, public SimulatorAccess<dim>
      {
        public:
          void
          apply_return_mapping(const MaterialModel::MaterialModelInputs<dim> &in,
                               const unsigned int                             i,
                               const std::vector<double> &                    volume_fractions,
                               const SymmetricTensor<2,dim> &                 trial_stress,
                               const SymmetricTensor<4,dim> &                 stiffness_tensor,
                               StateUpdates<dim> &                            state_updates) const override;

          void
          compute_tangent_modulus(const MaterialModel::MaterialModelInputs<dim> &in,
                                  const unsigned int                             i,
                                  const std::vector<double> &                    volume_fractions,
                                  const SymmetricTensor<2,dim> &                 trial_stress,
                                  const SymmetricTensor<4,dim> &                 stiffness_tensor,
                                  SymmetricTensor<4,dim> &                       tangent_modulus) const override;

          static
          void
          declare_parameters (ParameterHandler &prm);

          void
          parse_parameters (ParameterHandler &prm) override;

        private:
          std::pair<double, double>
          cohesion_and_hardening_modulus(const double plastic_strain,
                                         const double plastic_multiplier,
                                         const unsigned int field) const;

          std::vector<double> eta;
          std::vector<double> eta_bar;
          std::vector<double> xi;
          std::vector<double> c0;

          std::vector<double> start_weakening_intervals;
          std::vector<double> end_weakening_intervals;
          std::vector<double> weakening_factors;

          std::vector<double> damper_viscosity;
      };
    }
  }
}

#endif
