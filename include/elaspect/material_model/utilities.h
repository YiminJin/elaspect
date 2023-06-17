#ifndef _elaspect_material_model_utilities_h
#define _elaspect_material_model_utilities_h

#include <elaspect/global.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/table.h>

namespace elaspect
{
  namespace MaterialUtilities
  {
    using namespace dealii;

    std::vector<double>
    compute_composition_fractions(const std::vector<double> &compositional_fields);

    enum CompositionalAveragingOperation
    {
      harmonic,
      arithmetic,
      geometric,
      maximum_composition
    };

    CompositionalAveragingOperation
    parse_compositional_averaging_operation (const std::string &parameter_name,
                                             const ParameterHandler &prm);

    double average_value (const std::vector<double> &volume_fractions,
                          const std::vector<double> &parameter_values,
                          const CompositionalAveragingOperation &average_type);

    template <int dim>
    void 
    compute_isotropic_stiffness_tensor (const double            bulk_modulus,
                                        const double            shear_modulus,
                                        SymmetricTensor<4,dim> &stiffness_tensor);

    template <int dim>
    void
    compute_isotropic_compliance_tensor (const double            bulk_modulus,
                                         const double            shear_modulus,
                                         SymmetricTensor<4,dim> &compliance_tensor);

    template <int dim>
    double
    get_bulk_modulus(const SymmetricTensor<4,dim> &stiffness_tensor);

    template <int dim>
    double 
    get_shear_modulus(const SymmetricTensor<4,dim> &stiffness_tensor);

    template <int dim>
    SymmetricTensor<2,dim>
    self_outer_product (const Tensor<1,dim> &e);

    template <int dim>
    void
    compute_isotropic_tensor_derivative (const SymmetricTensor<2,dim> &X,
                                         const SymmetricTensor<2,dim> &y,
                                         const Table<2,double>        &dydx,
                                         SymmetricTensor<4,dim>       &D);
  }
}

#endif
