#ifndef _elaspect_material_model_io_interface_h
#define _elaspect_material_model_io_interface_h

#include <elaspect/global.h>
#include <elaspect/quadrature_point_data.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/numerics/data_postprocessor.h>

namespace elaspect
{
  template <int dim> struct Introspection;

  namespace QuadraturePointProperty
  {
    template <int dim> class Manager;
  }

  namespace MaterialModel
  {
    using namespace dealii;

    namespace FieldDependences
    {
      enum Dependence
      {
        none                    = 0,

        position                = 0x0001,
        displacement            = 0x0002,
        displacement_gradient   = 0x0004,
        fluid_pressure          = 0x0008,
        temperature             = 0x0010,
        composition             = 0x0020,
        old_stress              = 0x0040,

        all_dependences        = position |
                                 displacement |
                                 displacement_gradient |
                                 fluid_pressure |
                                 temperature |
                                 composition |
                                 old_stress
      };

      inline Dependence operator | (const Dependence d1,
                                    const Dependence d2)
      {
        return static_cast<Dependence>(static_cast<unsigned int>(d1) |
                                       static_cast<unsigned int>(d2));
      }

      inline Dependence &operator|=(Dependence &d1,
                                    const Dependence d2)
      {
        d1 = d1 | d2;
        return d1;
      }

      inline Dependence operator & (const Dependence d1,
                                    const Dependence d2)
      {
        return static_cast<Dependence>(static_cast<unsigned int>(d1) &
                                       static_cast<unsigned int>(d2));
      }

      inline Dependence &operator &= (Dependence &d1,
                                      const Dependence d2)
      {
        d1 = d1 & d2;
        return d1;
      }
    }

    namespace MaterialProperties
    {
      enum Property
      {
        none                            = 0,

        tangent_modulus                 = 0x0001,
        viscosity                       = 0x0002,
        density                         = 0x0004,
        specific_heat                   = 0x0008,
        thermal_expansivity             = 0x0010,
        thermal_conductivity            = 0x0020,
        entropy_derivative_temperature  = 0x0040,
        entropy_derivative_pressure     = 0x0080,
        reaction_terms                  = 0x0100,
        additional_outputs              = 0x0200,

        basic_properties                = tangent_modulus |
                                          density |
                                          specific_heat |
                                          thermal_expansivity |
                                          thermal_conductivity |
                                          reaction_terms,

        all_properties                  = viscosity |
                                          entropy_derivative_temperature |
                                          entropy_derivative_pressure |
                                          additional_outputs |
                                          basic_properties
      };

      inline Property operator | (const Property p1,
                                  const Property p2)
      {
        return static_cast<Property>(static_cast<unsigned int>(p1) |
                                     static_cast<unsigned int>(p2));
      }

      inline Property &operator|=(Property &p1,
                                  const Property p2)
      {
        p1 = p1 | p2;
        return p1;
      }

      inline Property operator & (const Property p1,
                                  const Property p2)
      {
        return static_cast<Property>(static_cast<unsigned int>(p1) &
                                     static_cast<unsigned int>(p2));
      }

      inline Property &operator&=(Property &p1,
                                  const Property p2)
      {
        p1 = p1 & p2;
        return p1;
      }
    }

    template <int dim>
    struct MaterialModelInputs
    {
      MaterialModelInputs(const unsigned int                  n_points,
                          const unsigned int                  n_comp,
                          const FieldDependences::Dependence dependences,
                          const MaterialProperties::Property  requested_properties);

      MaterialModelInputs(const DataPostprocessorInputs::Vector<dim> &input_data,
                          const Introspection<dim>                   &introspection,
                          const FieldDependences::Dependence         dependences,
                          const MaterialProperties::Property          requested_properties);

      void reinit(const FEValuesBase<dim>                   &fe_values,
                  const QPDHandler<dim>                     &qpd_handler,
                  const Introspection<dim>                  &introspection,
                  const TrilinosWrappers::MPI::BlockVector  &solution);

      MaterialModelInputs (const MaterialModelInputs &source) = default;

      unsigned int n_evaluation_points() const;

      std::vector<Point<dim>> position;

      std::vector<Tensor<1,dim>> incremental_displacement;

      std::vector<Tensor<2,dim>> incremental_displacement_gradient;

      std::vector<double> fluid_pressure;

      std::vector<double> temperature;

      std::vector<std::vector<double>> composition;

      std::vector<SymmetricTensor<2,dim>> old_stress;

      typename QPDHandler<dim>::active_cell_iterator qpd_cell;

      const FieldDependences::Dependence field_dependences;

      const MaterialProperties::Property requested_properties;
    };

    template <int dim> class AdditionalMaterialOutputs;

    template <int dim>
    struct MaterialModelOutputs
    {
      MaterialModelOutputs (const unsigned int n_points,
                            const unsigned int n_comp);

      MaterialModelOutputs (const MaterialModelOutputs &source);

      MaterialModelOutputs (MaterialModelOutputs &&) = default;

      MaterialModelOutputs &operator= (const MaterialModelOutputs &source) = delete;

      MaterialModelOutputs &operator= (MaterialModelOutputs &&) = default;

      std::vector<SymmetricTensor<4,dim> > tangent_moduli;

      std::vector<double> viscosities;

      std::vector<double> densities;

      std::vector<double> specific_heat;

      std::vector<double> thermal_expansivities;

      std::vector<double> thermal_conductivities;

      std::vector<double> entropy_derivative_temperature;

      std::vector<double> entropy_derivative_pressure;

      std::vector<std::vector<double> > reaction_terms;

      std::vector<std::unique_ptr<AdditionalMaterialOutputs<dim> > > additional_outputs;

      template <class AdditionalOutputType>
      AdditionalOutputType *get_additional_output();

      template <class AdditionalOutputType>
      const AdditionalOutputType *get_additional_output() const;
    };

    template <int dim>
    class AdditionalMaterialOutputs
    {
      public:
        virtual ~AdditionalMaterialOutputs () = default;
    };


// --------------------- template function definitions ----------------------------------

    template <int dim>
    template <class AdditionalOutputType>
    AdditionalOutputType *MaterialModelOutputs<dim>::get_additional_output()
    {
      for (unsigned int i=0; i<additional_outputs.size(); ++i)
        {
          AdditionalOutputType *result = dynamic_cast<AdditionalOutputType *> (additional_outputs[i].get());
          if (result)
            return result;
        }
      return nullptr;
    }


    template <int dim>
    template <class AdditionalOutputType>
    const AdditionalOutputType *MaterialModelOutputs<dim>::get_additional_output() const
    {
      for (unsigned int i=0; i<additional_outputs.size(); ++i)
        {
          const AdditionalOutputType *result = dynamic_cast<const AdditionalOutputType *> (additional_outputs[i].get());
          if (result)
            return result;
        }
      return nullptr;
    }
  }
}

#endif
