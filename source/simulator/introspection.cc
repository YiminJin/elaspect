#include <elaspect/introspection.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>

namespace elaspect
{
  namespace internal
  {
    template <int dim>
    typename Introspection<dim>::ComponentIndices
    setup_component_indices (const FEVariableCollection<dim> &fevs)
    {
      typename Introspection<dim>::ComponentIndices ci;
      for (unsigned int i = 0; i < dim; ++i)
        ci.displacement[i] = fevs.variable("displacement").first_component_index + i;

      ci.temperature = fevs.variable("temperature").first_component_index;

      ci.compositional_fields.clear();
      unsigned int n_compositional_fields = fevs.variable("composition").n_components();
      for (unsigned int i = 0; i < n_compositional_fields; ++i)
        ci.compositional_fields.push_back(fevs.variable("composition").first_component_index+i);

      ci.stress.clear();
      for (unsigned int i = 0; i < SymmetricTensor<2,dim>::n_independent_components; ++i)
        ci.stress.push_back(fevs.variable("stress").first_component_index+i);

      return ci;
    }


    template <int dim>
    typename Introspection<dim>::BlockIndices
    setup_blocks (const FEVariableCollection<dim> &fevs)
    {
      typename Introspection<dim>::BlockIndices bi;
      bi.displacement = fevs.variable("displacement").block_index;
      bi.temperature = fevs.variable("temperature").block_index;

      bi.compositional_fields.clear();
      unsigned int n_compositional_fields = fevs.variable("composition").n_components();
      for (unsigned int i = 0; i < n_compositional_fields; ++i)
        bi.compositional_fields.push_back(fevs.variable("composition").block_index+i);

      bi.stress.clear();
      for (unsigned int i = 0; i < SymmetricTensor<2,dim>::n_independent_components; ++i)
        bi.stress.push_back(fevs.variable("stress").block_index+i);

      return bi;
    }


    template <int dim>
    typename Introspection<dim>::BaseElements
    setup_base_elements (const FEVariableCollection<dim> &fevs)
    {
      typename Introspection<dim>::BaseElements base_elements;
      base_elements.displacement = fevs.variable("displacement").base_index;
      base_elements.temperature = fevs.variable("temperature").base_index;
      base_elements.compositional_fields = fevs.variable("composition").base_index;
      base_elements.stress = fevs.variable("stress").base_index;

      return base_elements;
    }


    template <int dim>
    unsigned int
    count_qpd_components(const Parameters<dim> &parameters)
    {
      unsigned int n_qpd_components = parameters.n_compositional_fields + // composition
                                      SymmetricTensor<2,dim>::n_independent_components; // old stress
      if (parameters.constitutive_relation & ConstitutiveRelation::plasticity)
        n_qpd_components += SymmetricTensor<2,dim>::n_independent_components // current stress
                            + 3; // old plastic strain, current plastic strain and return-mapping flag

      if (parameters.constitutive_relation & ConstitutiveRelation::thermal_expansion)
        n_qpd_components += 1; // old temperature

      return n_qpd_components;
    }


    template <int dim>
    std::shared_ptr<FiniteElement<dim>>
    new_FE_Q_or_DGQ(const bool discontinuous,
                    const unsigned int degree)
    {
      if (discontinuous)
        return std::make_shared<FE_DGQ<dim>>(degree);
      else
        return std::make_shared<FE_Q<dim>>(degree);
    }
  }


  template <int dim>
  std::vector<VariableDeclaration<dim> >
  construct_default_variables (const Parameters<dim> &parameters)
  {
    std::vector<VariableDeclaration<dim> > variables;

    variables.push_back(
      VariableDeclaration<dim>("displacement",
                               std::make_shared<FE_Q<dim>>(parameters.displacement_degree),
                               dim,
                               1));

    variables.push_back(
      VariableDeclaration<dim>("temperature",
                               std::make_shared<FE_Q<dim>>(parameters.temperature_degree),
                               1,
                               1));

    variables.push_back(
      VariableDeclaration<dim>("composition",
                               internal::new_FE_Q_or_DGQ<dim>(parameters.use_ALE_method,
                                                              parameters.composition_degree),
                               parameters.n_compositional_fields,
                               parameters.n_compositional_fields));

    const unsigned int displacement_gradient_degree = 
      (parameters.displacement_degree > 1 ? parameters.displacement_degree - 1 : 1);
    variables.push_back(
      VariableDeclaration<dim>("stress",
                               internal::new_FE_Q_or_DGQ<dim>(parameters.use_ALE_method,
                                                              displacement_gradient_degree),
                               SymmetricTensor<2,dim>::n_independent_components,
                               SymmetricTensor<2,dim>::n_independent_components));

    if (parameters.constitutive_relation & ConstitutiveRelation::plasticity)
      variables.push_back(
        VariableDeclaration<dim>("plastic strain",
                                 internal::new_FE_Q_or_DGQ<dim>(parameters.use_ALE_method,
                                                                displacement_gradient_degree),
                                 1,
                                 1));

    if (parameters.constitutive_relation & ConstitutiveRelation::pore_fluid)
      variables.push_back(
        VariableDeclaration<dim>("fluid pressure",
                                 std::make_shared<FE_Q<dim>>(parameters.fluid_pressure_degree),
                                 1,
                                 1));

    return variables;
  }


  template <int dim>
  Introspection<dim>::
  Introspection (const std::vector<VariableDeclaration<dim> > &variable_definition,
                 const Parameters<dim>                        &parameters)
    :
    FEVariableCollection<dim>(variable_definition),
    n_components(FEVariableCollection<dim>::n_components()),
    n_compositional_fields(parameters.n_compositional_fields),
    n_quadrature_points(Utilities::fixed_power<dim, int>(parameters.n_gaussian_points)),
    component_indices(internal::setup_component_indices<dim>(*this)),
    n_blocks(FEVariableCollection<dim>::n_blocks()),
    block_indices(internal::setup_blocks<dim>(*this)),
    extractors(component_indices),
    base_elements(internal::setup_base_elements<dim>(*this)),
    component_masks(*this),
    qpd_indicators(parameters),
    n_qpd_components(internal::count_qpd_components(parameters)),
    system_dofs_per_block(n_blocks),
    composition_names(parameters.names_of_compositional_fields)
  {}


  template <int dim>
  Introspection<dim>::~Introspection ()
  {}


  namespace
  {
    std::vector<FEValuesExtractors::Scalar>
    make_extractor_sequence(const std::vector<unsigned int> &component_indices)
    {
      std::vector<FEValuesExtractors::Scalar> x;
      for (const unsigned int component : component_indices)
        x.emplace_back(component);
      return x;
    }
  }


  template <int dim>
  Introspection<dim>::Extractors::
  Extractors (const Introspection<dim>::ComponentIndices &component_indices)
    :
    displacement (component_indices.displacement[0]),
    temperature (component_indices.temperature),
    compositional_fields (make_extractor_sequence(component_indices.compositional_fields)),
    stress (component_indices.stress[0])
  {}


  namespace
  {
    template <int dim>
    std::vector<ComponentMask>
    make_component_mask_sequence(const FEVariable<dim> &variable)
    {
      std::vector<ComponentMask> result;
      for (unsigned int i = 0; i < variable.multiplicity; ++i)
      {
        result.push_back(ComponentMask(variable.component_mask.size(), false));
        result.back().set(variable.first_component_index+i, true);
      }
      return result;
    }
  }


  template <int dim>
  Introspection<dim>::ComponentMasks::
  ComponentMasks (FEVariableCollection<dim> &fevs)
    :
    displacement (fevs.variable("displacement").component_mask),
    temperature (fevs.variable("temperature").component_mask),
    compositional_fields (make_component_mask_sequence(fevs.variable("composition"))),
    stress (make_component_mask_sequence(fevs.variable("stress")))
  {}


  template <int dim>
  Introspection<dim>::QPDIndicators::
  QPDIndicators (const Parameters<dim> &parameters)
  {
    composition = 0;
    old_stress  = parameters.n_compositional_fields;

    if (parameters.constitutive_relation & ConstitutiveRelation::plasticity)
    {
      current_stress         = old_stress + SymmetricTensor<2,dim>::n_independent_components;
      old_plastic_strain     = current_stress + SymmetricTensor<2,dim>::n_independent_components;
      current_plastic_strain = old_plastic_strain + 1;
      return_mapping_flag    = current_plastic_strain + 1;
    }
    else
    {
      current_stress         = numbers::invalid_unsigned_int;
      old_plastic_strain     = numbers::invalid_unsigned_int;
      current_plastic_strain = numbers::invalid_unsigned_int;
      return_mapping_flag    = numbers::invalid_unsigned_int;
    }
  }


  template <int dim>
  unsigned int
  Introspection<dim>::compositional_index_for_name(const std::string &name) const
  {
    const std::vector<std::string>::const_iterator
    it = std::find(composition_names.begin(), composition_names.end(), name);
    AssertThrow(it != composition_names.end(),
                ExcMessage("The compositional field " + name +
                           " you asked for is not used in the simulation."));

    return (it - composition_names.begin());
  }


  template <int dim>
  std::string
  Introspection<dim>::name_for_compositional_index(const unsigned int index) const
  {
    // make sure that what we get here is really an index of one of the compositional fields
    AssertIndexRange(index,composition_names.size());
    return composition_names[index];
  }
}


// explicit instantiations
namespace elaspect
{
#define INSTANTIATE(dim) \
  template struct Introspection<dim>; \
  template \
  std::vector<VariableDeclaration<dim> > \
  construct_default_variables (const Parameters<dim> &parameters);

  ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
