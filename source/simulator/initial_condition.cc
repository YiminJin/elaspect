#include <elaspect/simulator.h>
#include <elaspect/geometry_model/sphere.h>
#include <elaspect/geometry_model/spherical_shell.h>

#include <deal.II/numerics/vector_tools.h>

namespace elaspect
{
  template <int dim>
  void Simulator<dim>::initialize_temperature_field ()
  {
    // create a fully distributed vector since we
    // need to write into it and we can not
    // write into vectors with ghost elements
    TrilinosWrappers::MPI::BlockVector initial_solution;
    initial_solution.reinit (system_rhs, false);

    const unsigned int base_element = introspection.base_elements.temperature;
    const std::vector<Point<dim> > support_points
      = finite_element.base_element(base_element).get_unit_support_points();

    FEValues<dim> fe_values (*mapping, finite_element,
                             support_points,
                             update_quadrature_points);

    std::vector<types::global_dof_index> local_dof_indices (finite_element.dofs_per_cell);

    const VectorFunctionFromScalarFunctionObject<dim, double> T_init_function 
      (
        [&](const Point<dim> &p) -> double
        {
          return initial_temperature_manager.initial_temperature(p);
        },
        introspection.component_indices.temperature,
        introspection.n_components
      );

    VectorTools::interpolate (*mapping,
                              dof_handler,
                              T_init_function,
                              initial_solution,
                              introspection.component_masks.temperature);

    initial_solution.compress (VectorOperation::insert);

    // then apply constraints and copy the
    // result into vectors with ghost elements. to do so,
    // we need the current constraints to be correct for
    // the current time
    compute_current_constraints();
    current_constraints.distribute(initial_solution);

    // Now copy the temperature block into the solution variables
    const unsigned int block_idx = introspection.block_indices.temperature;

    solution.block(block_idx) = initial_solution.block(block_idx);
    old_solution.block(block_idx) = initial_solution.block(block_idx);
  }


  template <int dim>
  void Simulator<dim>::initialize_fluid_pressure_field ()
  {

  }


  template <int dim>
  void Simulator<dim>::initialize_compositional_fields()
  {
    // In this function, we shall initialize not only the quadrature point data,
    // but also the solution vector blocks corresponding to compositional fields,
    // since some mesh refinement strategies (like composition gradient) may 
    // rely on the solution vector.

    // Initialize the quadrature point data.
    const unsigned int n_q_points = qpd_handler.n_quadrature_points();
    std::vector<Point<dim>> quadrature_points(n_q_points);
    std::vector<double> composition_values(n_q_points);

    for (const auto &cell : triangulation.active_cell_iterators())
      if (cell->is_locally_owned())
      {
        const typename QPDHandler<dim>::active_cell_iterator qpd_cell(*cell, &qpd_handler);

        for (unsigned int c = 0; c < introspection.n_compositional_fields; ++c)
        {
          for (unsigned int q = 0; q < n_q_points; ++q)
          {
            const Point<dim> position = mapping->transform_unit_to_real_cell(cell, quadrature_formula.point(q));
            composition_values[q] = initial_composition_manager.initial_composition(position, c);
          }

          qpd_cell->set(introspection.qpd_indicators.composition + c, composition_values);
        }
      }

    // Initialize the solution vector.
    TrilinosWrappers::MPI::BlockVector initial_solution(introspection.index_sets.system_partitioning, 
                                                        mpi_communicator);

    for (unsigned int c = 0; c < introspection.n_compositional_fields; ++c)
    {
      const VectorFunctionFromScalarFunctionObject<dim> comp_init_function(
        [&] (const Point<dim> &x) -> double
        {
          return initial_composition_manager.initial_composition(x, c);
        },
        introspection.component_indices.compositional_fields[c],
        introspection.n_components);

      VectorTools::interpolate(*mapping,
                               dof_handler,
                               comp_init_function,
                               initial_solution,
                               introspection.component_masks.compositional_fields[c]);

      initial_solution.compress(VectorOperation::insert);
    }

    compute_current_constraints();
    current_constraints.distribute(initial_solution);

    for (unsigned int c = 0; c < introspection.n_compositional_fields; ++c)
    {
      const unsigned int block_idx = introspection.block_indices.compositional_fields[c];
      solution.block(block_idx) = initial_solution.block(block_idx);
      old_solution.block(block_idx) = initial_solution.block(block_idx);
    }
  }
}


// explicit instantiations
namespace elaspect
{
#define INSTANTIATE(dim) \
  template void Simulator<dim>::initialize_temperature_field(); \
  template void Simulator<dim>::initialize_fluid_pressure_field(); \
  template void Simulator<dim>::initialize_compositional_fields();

  ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE 
}
