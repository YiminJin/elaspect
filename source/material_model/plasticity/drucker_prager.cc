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


#include <elaspect/material_model/plasticity/drucker_prager.h>
#include <elaspect/utilities.h>

#include <deal.II/physics/elasticity/standard_tensors.h>

namespace elaspect
{
  namespace MaterialModel
  {
    namespace Plasticity
    {
      using namespace dealii::Physics::Elasticity;

      template <int dim>
      void
      DruckerPrager<dim>::
      apply_return_mapping(const MaterialModel::MaterialModelInputs<dim> &in,
                           const unsigned int                             i,
                           const std::vector<double> &                    volume_fractions,
                           const SymmetricTensor<2,dim> &                 trial_stress,
                           const SymmetricTensor<4,dim> &                 stiffness_tensor,
                           StateUpdates<dim> &                            state_updates) const
      {
        const SymmetricTensor<2,dim> s_tr = deviator(trial_stress);
        const double sqrt_J2_tr = std::sqrt(std::abs(second_invariant(s_tr)));
        const double p_tr = trace(trial_stress) / dim;

        const double epsilon_p_old = std::max(
          in.qpd_cell->get_scalar(i, this->introspection().qpd_indicators.old_plastic_strain), 0.0);

        // check plastic admissibility
        const unsigned int field = std::distance(volume_fractions.begin(),
                                                 std::max_element(volume_fractions.begin(),
                                                                  volume_fractions.end()));

        std::pair<double, double> c_and_H = 
          cohesion_and_hardening_modulus(epsilon_p_old, 0, field);
        double F = sqrt_J2_tr + eta[field] * p_tr - xi[field] * c_and_H.first;
        if (F <= 0)
        {
          state_updates.stress         = trial_stress;
          state_updates.flag           = 0;
          state_updates.plastic_strain = epsilon_p_old;
          return;
        }

        const double G = MaterialUtilities::get_shear_modulus(stiffness_tensor),
                     K = MaterialUtilities::get_bulk_modulus(stiffness_tensor);

        const double tol = 1.e-8;
        const unsigned int iter_max = 10;

        double epsilon_p = numbers::signaling_nan<double>();

        // set initial guess for dgamma
        double dgamma = 0;

        // apply return mapping to the smooth wall
        for (unsigned int iter = 0; iter < iter_max; ++iter)
        {
          double d = -(G + K * eta[field] * eta_bar[field] +
                       xi[field] * xi[field] * c_and_H.second);
          AssertThrow(d < 0, ExcInternalError());
          dgamma -= F / d;

          // compute new residual
          epsilon_p = epsilon_p_old + xi[field] * dgamma;
          c_and_H = cohesion_and_hardening_modulus(epsilon_p, dgamma, field);
          F = sqrt_J2_tr - G * dgamma + 
              eta[field] * (p_tr - K * eta_bar[field] * dgamma) - 
              xi[field] * c_and_H.first;
          if (std::abs(F) < tol * c_and_H.first)
            break;
        }

        // check validity
        if (sqrt_J2_tr - G * dgamma >= 0)
        {
          // smooth wall return is valid -- update stress and other variables
          state_updates.stress = (1 - G * dgamma / sqrt_J2_tr) * s_tr;
          const double p_new = p_tr - K * eta_bar[field] * dgamma;
          for (unsigned int d = 0; d < dim; ++d)
            state_updates.stress[d][d] += p_new;

          // for convenience, we store the value of plastic multiplier in 
          // QPD slot 'flag' (since the plastic multiplier is always positive), 
          // so that we do not have to calculate the plastic multiplier 
          // when calculating the tangent modulus
          state_updates.flag = dgamma;
          state_updates.plastic_strain = epsilon_p;
        }
        else
        {
          // smooth wall return not valid -- return to apex
          const double alpha = xi[field] / eta_bar[field];
          const double beta  = xi[field] / eta[field];

          // set initial guess for depsilon_vp
          double depsilon_vp = 0;

          c_and_H = cohesion_and_hardening_modulus(epsilon_p_old, 0, field);
          double res = beta * c_and_H.first - p_tr;

          for (unsigned int iter = 0; iter < iter_max; ++iter)
          {
            double d = alpha * beta * c_and_H.second + K;
            AssertThrow(d > 0, ExcInternalError());
            depsilon_vp -= res / d;

            // compute new residual
            epsilon_p = epsilon_p_old + alpha * depsilon_vp;
            c_and_H = cohesion_and_hardening_modulus(epsilon_p, 
                                                     depsilon_vp / eta_bar[field], 
                                                     field);
            res = beta * c_and_H.first - (p_tr - K * depsilon_vp);
            if (std::abs(res) < tol * c_and_H.first)
              break;
          }

          const double p_new = p_tr - K * depsilon_vp;
          state_updates.stress = 0;
          for (unsigned int d = 0; d < dim; ++d)
            state_updates.stress[d][d] = p_new;

          state_updates.flag = -1;
          state_updates.plastic_strain = epsilon_p;
        }
      }


      template <int dim>
      void
      DruckerPrager<dim>::
      compute_tangent_modulus(const MaterialModel::MaterialModelInputs<dim> &in,
                              const unsigned int                             i,
                              const std::vector<double> &                    volume_fractions,
                              const SymmetricTensor<2,dim> &                 trial_stress,
                              const SymmetricTensor<4,dim> &                 stiffness_tensor,
                              SymmetricTensor<4,dim> &                       tangent_modulus) const
      {
        const double flag = in.qpd_cell->get_scalar(
          i, this->introspection().qpd_indicators.return_mapping_flag);
        const double epsilon_p = in.qpd_cell->get_scalar(
          i, this->introspection().qpd_indicators.current_plastic_strain);

        if (flag == 0)
        {
          // stress is inside the Drucker-Prager cone. the tangent modulus is
          // equal to the stiffness modulus.
          tangent_modulus = stiffness_tensor;
          return;
        }

        // stress is outside the Drucker-Prager cone. calculate the consistent 
        // tangent modulus.
        const double G = MaterialUtilities::get_shear_modulus(stiffness_tensor),
                     K = MaterialUtilities::get_bulk_modulus(stiffness_tensor);

        const unsigned int field = std::distance(volume_fractions.begin(),
                                                 std::max_element(volume_fractions.begin(),
                                                                  volume_fractions.end()));

        // here we do not need the value of cohesion, so the second argument
        // of function cohesion_and_hardening_modulus does not affect the 
        // final result
        std::pair<double, double> c_and_H = cohesion_and_hardening_modulus(epsilon_p, 0, field);

        if (flag == -1)
        {
          // the stress has been mapped to the apex.
          const double alpha = xi[field] / eta_bar[field];
          const double beta  = xi[field] / eta[field];

          tangent_modulus = (K * (1 - K / (K + alpha * beta * c_and_H.second)))
                            * StandardTensors<dim>::IxI;
        }
        else
        {
          // the stress has been mapped to the smooth portion of the cone.
          const double dgamma = flag;
          AssertThrow(dgamma > 0,
                      ExcMessage("The plastic multiplier is non-positive."));

          const SymmetricTensor<2,dim> s_tr = deviator(trial_stress);
          const double sqrt_J2_tr = std::sqrt(std::fabs(second_invariant(s_tr)));
          const SymmetricTensor<2,dim> N = s_tr * (numbers::SQRT1_2 / sqrt_J2_tr);

          const double A_inv = G + K * eta[field] * eta_bar[field] + 
                               xi[field] * xi[field] * c_and_H.second;
          AssertThrow(A_inv > 0,
                      ExcMessage("The absolute value of softening modulus "
                                 "($-dc/d\\bar{\\varepsilon}^p$) is too large."));

          const double A = 1. / A_inv;
          const double B = dgamma / sqrt_J2_tr;

          tangent_modulus = 2 * G * 
                            ( ( 1 - B * G ) * StandardTensors<dim>::dev_P +
                              ( ( B - A ) * G) * outer_product(N, N) )
                            -
                            numbers::SQRT2 * G * A * K *
                            ( outer_product(N, StandardTensors<dim>::I) * eta[field] +
                              outer_product(StandardTensors<dim>::I, N) * eta_bar[field] )
                            +
                            ( K * ( 1 - K * eta[field] * eta_bar[field] * A ) ) *
                            StandardTensors<dim>::IxI;
        }
      }


      template <int dim>
      std::pair<double, double>
      DruckerPrager<dim>::
      cohesion_and_hardening_modulus(const double plastic_strain,
                                     const double plastic_multiplier,
                                     const unsigned int field) const
      {
        // viscous regularization
        const double R = damper_viscosity[field] / (this->get_timestep() * xi[field]);

        if (plastic_strain < start_weakening_intervals[field])
          return std::make_pair(c0[field] + R * plastic_multiplier, R);

        if (plastic_strain > end_weakening_intervals[field])
          return std::make_pair(c0[field] * weakening_factors[field] + R * plastic_multiplier, R);

        const double H = c0[field] * (weakening_factors[field] - 1.0) /
                         (end_weakening_intervals[field] - start_weakening_intervals[field]);

        const double c = c0[field] + H * (plastic_strain - start_weakening_intervals[field]);

        return std::make_pair(c + R * plastic_multiplier, H + R);
      }


      template <int dim>
      void
      DruckerPrager<dim>::declare_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Material model");
        {
          prm.enter_subsection("Plasticity");
          {
            prm.enter_subsection("Drucker Prager");
            {
              prm.declare_entry("Friction angle", "30",
                                Patterns::List(Patterns::Double(0, 45)),
                                "");
              prm.declare_entry("Dilatancy angle", "30",
                                Patterns::List(Patterns::Double(0, 45)),
                                "");
              prm.declare_entry("Cohesion", "1e6",
                                Patterns::List(Patterns::Double(0)),
                                "");
              prm.declare_entry("Start weakening intervals", "0",
                                Patterns::List(Patterns::Double(0)),
                                "");
              prm.declare_entry("End weakening intervals", "1",
                                Patterns::List(Patterns::Double(0)),
                                "");
              prm.declare_entry("Weakening factors", "1",
                                Patterns::List(Patterns::Double(0)),
                                "");
              prm.declare_entry("Plastic damper viscosity", "0",
                                Patterns::List(Patterns::Double(0)),
                                "");
              prm.declare_entry("Use compression cone", "true",
                                Patterns::Bool(),
                                "");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      DruckerPrager<dim>::parse_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Material model");
        {
          prm.enter_subsection("Plasticity");
          {
            prm.enter_subsection("Drucker Prager");
            {
              const unsigned int n_fields = this->n_compositional_fields() + 1;

              std::vector<double> phi, psi;
              phi = Utilities::possibly_extend_from_1_to_N(
                Utilities::string_to_double(Utilities::split_string_list(prm.get("Friction angle"))),
                n_fields,
                "Friction angle");
              psi = Utilities::possibly_extend_from_1_to_N(
                Utilities::string_to_double(Utilities::split_string_list(prm.get("Dilatancy angle"))),
                n_fields,
                "Dilatancy angle");
              c0 = Utilities::possibly_extend_from_1_to_N(
                Utilities::string_to_double(Utilities::split_string_list(prm.get("Cohesion"))),
                n_fields,
                "Cohesion");
              start_weakening_intervals = Utilities::possibly_extend_from_1_to_N(
                Utilities::string_to_double(Utilities::split_string_list(prm.get("Start weakening intervals"))),
                n_fields,
                "Start weakening intervals");
              end_weakening_intervals = Utilities::possibly_extend_from_1_to_N(
                Utilities::string_to_double(Utilities::split_string_list(prm.get("End weakening intervals"))),
                n_fields,
                "End weakening intervals");
              weakening_factors = Utilities::possibly_extend_from_1_to_N(
                Utilities::string_to_double(Utilities::split_string_list(prm.get("Weakening factors"))),
                n_fields,
                "Weakening factors");
              damper_viscosity = Utilities::possibly_extend_from_1_to_N(
                Utilities::string_to_double(Utilities::split_string_list(prm.get("Plastic damper viscosity"))),
                n_fields,
                "Plastic damper viscosity");

              const bool use_compression_cone = prm.get_bool("Use compression cone");

              eta.resize(n_fields);
              eta_bar.resize(n_fields);
              xi.resize(n_fields);

              for (unsigned int field = 0; field < n_fields; ++field)
              {
                AssertThrow(phi[field] > 0, 
                            ExcMessage("The Drucker Prager model requires the "
                                       "friction angle to be greater than zero."));
                AssertThrow(psi[field] > 0, 
                            ExcMessage("The Drucker Prager model requires the "
                                       "dilatancy angle to be greater than zero."));
                AssertThrow(start_weakening_intervals[field] <= end_weakening_intervals[field],
                            ExcMessage("Invalid weakening interval detected."));

                phi[field] *= numbers::PI / 180.;
                psi[field] *= numbers::PI / 180.;

                if (dim == 2)
                {
                  eta[field]     = std::sin(phi[field]);
                  eta_bar[field] = std::sin(psi[field]);
                  xi[field]      = std::cos(phi[field]);
                }
                else
                {
                  double factor1, factor2;
                  if (use_compression_cone)
                  {
                    factor1 = 1. / (std::sqrt(3.) * (3 - std::sin(phi[field])));
                    factor2 = 1. / (std::sqrt(3.) * (3 - std::sin(psi[field])));
                  }
                  else
                  {
                    factor1 = 1. / (std::sqrt(3.) * (3 + std::sin(phi[field])));
                    factor2 = 1. / (std::sqrt(3.) * (3 + std::sin(psi[field])));
                  }

                  eta[field]     = (6 * std::sin(phi[field])) * factor1; 
                  eta_bar[field] = (6 * std::sin(psi[field])) * factor2; 
                  xi[field]      = (6 * std::cos(phi[field])) * factor1;
                }
              }
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
  namespace MaterialModel
  {
    namespace Plasticity
    {
      ELASPECT_REGISTER_PLASTICITY_MODEL(DruckerPrager,
                                     "drucker prager",
                                     "")
    }
  }
}
