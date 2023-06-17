#include <elaspect/material_model/utilities.h>

#include <deal.II/physics/elasticity/standard_tensors.h>

namespace elaspect
{
  namespace MaterialUtilities
  {
    using namespace Physics::Elasticity;

    std::vector<double>
    compute_composition_fractions(const std::vector<double> &compositional_fields)
    {
      std::vector<double> composition_fractions(compositional_fields.size() + 1);

      // Clip the compositional fields so they are between zero and one,
      // and sum the compositional fields for normalization purposes.
      double sum_composition = 0.0;
      std::vector<double> x_comp = compositional_fields;
      for (unsigned int i = 0; i < x_comp.size(); ++i)
      {
        x_comp[i] = std::min(std::max(x_comp[i], 0.0), 1.0);
        sum_composition += x_comp[i];
      }

      // Compute background field fractions
      if (sum_composition >= 1.0)
        composition_fractions[0] = 0.0;
      else
        composition_fractions[0] = 1.0 - sum_composition;

      // Compute and possibly normalize field fractions
      for (unsigned int i = 0; i < x_comp.size(); ++i)
      {
        if (sum_composition >= 1.0)
          composition_fractions[i+1] = x_comp[i] / sum_composition;
        else
          composition_fractions[i+1] = x_comp[i];
      }

      return composition_fractions;
    }


    CompositionalAveragingOperation
    parse_compositional_averaging_operation (const std::string &parameter_name,
                                             const ParameterHandler &prm)
    {
      CompositionalAveragingOperation averaging_operation;
      if (prm.get (parameter_name) == "harmonic")
        averaging_operation = MaterialUtilities::harmonic;
      else if (prm.get (parameter_name) == "arithmetic")
        averaging_operation = MaterialUtilities::arithmetic;
      else if (prm.get (parameter_name) == "geometric")
        averaging_operation = MaterialUtilities::geometric;
      else if (prm.get (parameter_name) == "maximum composition")
        averaging_operation = MaterialUtilities::maximum_composition;
      else
        {
          AssertThrow(false, ExcMessage("Not a valid viscosity averaging scheme"));

          //We will never get here, but we have to return something so the compiler does not complain
          return MaterialUtilities::harmonic;
        }

      return averaging_operation;
    }


    double
    average_value (const std::vector<double> &volume_fractions,
                   const std::vector<double> &parameter_values,
                   const enum CompositionalAveragingOperation &average_type)
    {
      Assert(volume_fractions.size() == parameter_values.size(),
             ExcMessage ("The volume fractions and parameter values vectors used for averaging "
                         "have to have the same length!"));

      double averaged_parameter = 0.0;

      switch (average_type)
      {
        case arithmetic:
        {
          for (unsigned int i = 0; i < volume_fractions.size(); ++i)
            averaged_parameter += volume_fractions[i] * parameter_values[i];
          break;
        }
        case harmonic:
        {
          for (unsigned int i = 0; i < volume_fractions.size(); ++i)
          {
            AssertThrow(parameter_values[i] > 0,
                        ExcMessage ("All parameter values must be greater than 0 for harmonic averaging!"));
            averaged_parameter += volume_fractions[i] / (parameter_values[i]);
          }
          averaged_parameter = 1.0 / averaged_parameter;
          break;
        }
        case geometric:
        {
          for (unsigned int i = 0; i < volume_fractions.size(); ++i)
          {
            AssertThrow(parameter_values[i] > 0,
                        ExcMessage ("All parameter values must be greater than 0 for geometric averaging!"));
            averaged_parameter += volume_fractions[i] * std::log(parameter_values[i]);
          }
          averaged_parameter = std::exp(averaged_parameter);
          break;
        }
        case maximum_composition:
        {
          const unsigned int idx = static_cast<unsigned int>(std::max_element( volume_fractions.begin(),
                                                                               volume_fractions.end() )
                                                             - volume_fractions.begin());
          averaged_parameter = parameter_values[idx];
          break;
        }
        default:
        {
          AssertThrow(false, ExcNotImplemented());
          break;
        }
      }
      return averaged_parameter;
    }


    template <int dim>
    void
    compute_isotropic_stiffness_tensor (const double            K,
                                        const double            G,
                                        SymmetricTensor<4,dim> &D)
    {
      const double lambda = K - 2. / dim * G;

      D = 0;
      for (unsigned int i = 0; i < dim; ++i)
        for (unsigned int j = 0; j < dim; ++j)
          for (unsigned int k = 0; k < dim; ++k)
            for (unsigned int l = 0; l < dim; ++l)
              D[i][j][k][l] = ( ( (i == k) && (j == l) ? G : 0.0 ) +
                                ( (i == l) && (j == k) ? G : 0.0 ) +
                                ( (i == j) && (k == l) ? lambda : 0.0 ) );
    }


    template <int dim>
    void
    compute_isotropic_compliance_tensor (const double            K,
                                         const double            G,
                                         SymmetricTensor<4,dim> &C)
    {
      const double a = 1. / (4 * G);
      const double b = ( dim == 2 ?
                         (G - K) / (4 * K * G) :
                         //(2 * G - 3 * K) / (4 * G * (G + 3 * K)) :
                         (2 * G - 3 * K) / (18 * K * G) );

      C = 0;
      for (unsigned int i = 0; i < dim; ++i)
        for (unsigned int j = 0; j < dim; ++j)
          for (unsigned int k = 0; k < dim; ++k)
            for (unsigned int l = 0; l < dim; ++l)
              C[i][j][k][l] = ( ( (i == k) && (j == l) ? a : 0.0 ) +
                                ( (i == l) && (j == k) ? a : 0.0 ) +
                                ( (i == j) && (k == l) ? b : 0.0 ) );
    }


    template <int dim>
    double
    get_bulk_modulus(const SymmetricTensor<4,dim> &stiffness_tensor)
    {
      return stiffness_tensor[0][0][1][1] + stiffness_tensor[0][1][0][1] * 2. / dim;
    }


    template <int dim>
    double
    get_shear_modulus(const SymmetricTensor<4,dim> &stiffness_tensor)
    {
      return stiffness_tensor[0][1][0][1];
    }


    template <int dim>
    SymmetricTensor<2,dim>
    self_outer_product (const Tensor<1,dim> &e)
    {
      SymmetricTensor<2,dim> E;
      for (unsigned int i = 0; i < dim; ++i)
        for (unsigned int j = i; j < dim; ++j)
          E[i][j] = e[i] * e[j];

      return E;
    }

    namespace internal
    {
      void
      assemble_derivatives(const SymmetricTensor<2, 2>                &/*X*/,
                           const SymmetricTensor<2, 2>                &/*Y*/,
                           const std::array<double, 2>                &x,
                           const std::array<double, 2>                &y,
                           const std::array<SymmetricTensor<2, 2>, 2> &E,
                           const Table<2, double>                     &dydx,
                           SymmetricTensor<4,2>                       &D)
      {
        D = 0;

        //TODO: what if max(x0, x1) = 0?
        const double tol = std::max(std::abs(x[0]),
                                    std::abs(x[1])) * 1e-8;
        if (std::abs(x[0] - x[1]) > tol)
        {
          SymmetricTensor<4, 2> ExE[2][2];
          for (unsigned int i = 0; i < 2; ++i)
            for (unsigned int j = 0; j < 2; ++j)
              ExE[i][j] = outer_product(E[i], E[j]);

          const double factor = (y[0] - y[1]) / (x[0] - x[1]);
          D = (StandardTensors<2>::S - ExE[0][0] - ExE[1][1]) * factor;

          for (unsigned int i = 0; i < 2; ++i)
            for (unsigned int j = 0; j < 2; ++j)
              D += dydx[i][j] * ExE[i][j];
        }
        else
        {
          D = StandardTensors<2>::S * (dydx[0][0] - dydx[0][1]) +
              StandardTensors<2>::IxI * dydx[0][1];
        }
      }


      void
      assemble_derivatives(const SymmetricTensor<2, 3>                &X,
                           const SymmetricTensor<2, 3>                &/*Y*/,
                           const std::array<double, 3>                &x,
                           const std::array<double, 3>                &y,
                           const std::array<SymmetricTensor<2, 3>, 3> &E,
                           const Table<2, double>                     &dydx,
                           SymmetricTensor<4,3>                       &D)
      {
        D = 0;

        // compute dX^2 / dX
        SymmetricTensor<4, 3> dX2dX;
        for (unsigned int i = 0; i < 3; ++i)
          for (unsigned int j = 0; j < 3; ++j)
            for (unsigned int k = 0; k < 3; ++k)
              for (unsigned int l = 0; l < 3; ++l)
                dX2dX[i][j][k][l] = 0.5 * ( (i == k ? X[j][l] : 0) +
                                            (i == l ? X[j][k] : 0) +
                                            (j == l ? X[i][k] : 0) +
                                            (j == k ? X[i][l] : 0) );

        const double tol = std::max(std::max(std::abs(x[0]),
                                             std::abs(x[1])),
                                    std::abs(x[2])) * 1e-8;
        if (std::abs(x[0] - x[1]) > tol || std::abs(x[1] - x[2]) > tol)
        {
          // x0 != x1 != x2
          SymmetricTensor<4, 3> ExE[3][3];
          for (unsigned int i = 0; i < 3; ++i)
            for (unsigned int j = 0; j < 3; ++j)
              ExE[i][j] = outer_product(E[i], E[j]);

          for (unsigned int i = 0; i < 3; ++i)
          {
            const unsigned int a = i;
            const unsigned int b = (i+1)%3;
            const unsigned int c = (i+2)%3;

            const double factor = y[a] / ((x[a] - x[b]) * (x[a] - x[c]));
            D += factor * (dX2dX - (x[b] + x[c]) * StandardTensors<3>::S
                           - (x[a] - x[b] + x[a] - x[c]) * ExE[a][a]
                           - (x[b] - x[c]) * (ExE[b][b] - ExE[c][c]));
          }

          for (unsigned int i = 0; i < 3; ++i)
            for (unsigned int j = 0; j < 3; ++j)
              D += dydx[i][j] * ExE[i][j];

          return;
        }
        
        if (std::abs(x[0] - x[1]) > tol || std::abs(x[1] - x[2]) > tol)
        {
          unsigned int a, b, c;
          if (std::abs(x[0] - x[1]) <= tol)
          {
            // x0 == x1 !== x2
            a = 2;
            b = 0;
            c = 1;
          }
          else
          {
            // x0 != x1 == x2
            a = 0;
            b = 1;
            c = 2;
          }

          const double factor = 1. / (x[a] - x[c]);

          const double s1 = factor * (factor * (y[a] - y[c]) + (dydx[c][b] - dydx[c][c]));
          const double s2 = factor * (2 * factor * x[c] * (y[a] - y[c]) + (x[a] + x[c]) * (dydx[c][b] - dydx[c][c]));
          const double s3 = factor * factor * (2 * factor * (y[a] - y[c]) + (dydx[a][c] + dydx[c][a] - dydx[a][a] - dydx[c][c]));
          const double s4 = factor * (2 * factor * factor * x[c] * (y[a] - y[c]) + (dydx[a][c] - dydx[c][b]) + factor * x[c] * (dydx[a][c] + dydx[c][a] - dydx[a][a] - dydx[c][c]));
          const double s5 = s4 + factor * (dydx[c][a] - dydx[a][c]);
          const double s6 = factor * (2 * factor * factor * x[c] * x[c] * (y[a] - y[c]) + factor * x[a] * x[c] * (dydx[a][c] + dydx[c][a]) - factor * x[c] * x[c] * (dydx[a][a] + dydx[c][c]) - (x[a] + x[c]) * dydx[c][b]);

          const SymmetricTensor<2,3> I = unit_symmetric_tensor<3>();
          SymmetricTensor<2,3> X;
          for (unsigned int i = 0; i < 3; ++i)
            X[i][i] = x[i];

          D = s1 * dX2dX -
              s2 * identity_tensor<3>() -
              s3 * outer_product(X, X) +
              s4 * outer_product(X, I) +
              s5 * outer_product(I, X) -
              s6 * outer_product(I, I);

          return;
        }

        // x0 == x1 == x2
        D = (dydx[1][1] - dydx[1][2]) * StandardTensors<3>::S + dydx[1][2] * StandardTensors<3>::IxI;
      }
    }


    template <int dim>
    void
    compute_isotropic_tensor_derivative(const SymmetricTensor<2,dim> &X,
                                        const SymmetricTensor<2,dim> &Y,
                                        const Table<2,double>        &dydx,
                                        SymmetricTensor<4,dim>       &D)
    {
      // compute eigenvalues and eigenprojections of X
      const std::array<std::pair<double, Tensor<1, dim>>,
                       std::integral_constant<int, dim>::value>
      x_and_e = eigenvectors(X, SymmetricTensorEigenvectorMethod::hybrid);
      std::array<double, dim> x;
      std::array<SymmetricTensor<2, dim>, dim> E;
      for (unsigned int i = 0; i < dim; ++i)
      {
        x[i] = x_and_e[i].first;
        E[i] = self_outer_product(x_and_e[i].second);
      }

      // compute eigenvalues of Y
      const std::array<double, dim> y = eigenvalues(Y);

      // assemble the derivatives
      internal::assemble_derivatives(X, Y, x, y, E, dydx, D);
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace MaterialUtilities
  {
#define INSTANTIATE(dim) \
    template void \
    compute_isotropic_stiffness_tensor<dim> (const double, \
                                             const double, \
                                             SymmetricTensor<4,dim> &); \
    \
    template void \
    compute_isotropic_compliance_tensor<dim> (const double, \
                                              const double, \
                                              SymmetricTensor<4,dim> &); \
    \
    template double \
    get_bulk_modulus(const SymmetricTensor<4,dim> &stiffness_tensor); \
    \
    template double \
    get_shear_modulus(const SymmetricTensor<4,dim> &stiffness_tensor); \
    \
    template SymmetricTensor<2,dim> \
    self_outer_product<dim> (const Tensor<1,dim> &); \
    \
    template void \
    compute_isotropic_tensor_derivative<dim> (const SymmetricTensor<2,dim> &, \
                                              const SymmetricTensor<2,dim> &, \
                                              const Table<2,double> &, \
                                              SymmetricTensor<4,dim> &);

    ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
