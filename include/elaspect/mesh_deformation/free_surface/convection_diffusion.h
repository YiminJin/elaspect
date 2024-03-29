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


#ifndef _elaspect_mesh_deformation_free_surface_convection_diffusion_h
#define _elaspect_mesh_deformation_free_surface_convection_diffusion_h

#include <elaspect/mesh_deformation/free_surface/interface.h>
#include <elaspect/simulator_access.h>

#include <deal.II/base/function.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1_eulerian.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>

namespace elaspect
{
  namespace MeshDeformation
  {
    namespace FreeSurface
    {
      template <int dim>
      class ConvectionDiffusion : public Interface<dim>, public SimulatorAccess<dim>
      {
        public:
          ConvectionDiffusion();

          void update() override;

          void initialize() override;

          void
          make_boundary_constraints(const TrilinosWrappers::MPI::Vector &mesh_displacements,
                                    const DoFHandler<dim>               &mesh_deformation_dof_handler,
                                    AffineConstraints<double>           &mesh_deformation_constraints,
                                    const std::set<types::boundary_id>  &fs_boundary_ids) const override;

          static
          void 
          declare_parameters(ParameterHandler &prm);

          void
          parse_parameters(ParameterHandler &prm) override;

        private:
          static const unsigned int surface_dim = dim - 1;

          void create_interface_mesh();

          void gather_data_to_root_process();

          void project_velocities_onto_surface();

          void compute_topographic_evolution();

          void execute_mesh_smoothing();

          void interpolate_coordinates_onto_interface();

          void scatter_data_to_surface_owners();

          unsigned int subdivide_time_step() const;

          void update_topography_in_material_coordinates();

          double 
          compute_artificial_viscosity(const std::vector<double> &               old_z_values,
                                       const std::vector<double> &               old_old_z_values,
                                       const std::vector<Tensor<1,surface_dim>> &old_z_grads,
                                       const double                              old_z_laplacian,
                                       const std::vector<Tensor<1,surface_dim>> &velocity_values,
                                       const std::vector<double> &               uplift_rates,
                                       const double                              global_u_infty,
                                       const double                              global_z_variation,
                                       const double                              global_z_average,
                                       const double                              global_entropy_variation,
                                       const double                              step_size,
                                       const double                              cell_diameter) const;

          struct BoundaryFaceId
          {
            CellId       adjacent_cell;
            unsigned int face_no;
            unsigned int surface_mpi_rank;
   
            BoundaryFaceId ()
              : face_no(numbers::invalid_unsigned_int)
              , surface_mpi_rank(numbers::invalid_unsigned_int)
            {}

            BoundaryFaceId (const CellId &adjacent_cell_,
                            const unsigned int face_no_)
              : adjacent_cell(adjacent_cell_)
              , face_no(face_no_)
              , surface_mpi_rank(numbers::invalid_unsigned_int)
            {}

            BoundaryFaceId (const CellId &     adjacent_cell_,
                            const unsigned int face_no_,
                            const unsigned int surface_mpi_rank_)
              : adjacent_cell(adjacent_cell_)
              , face_no(face_no_)
              , surface_mpi_rank(surface_mpi_rank_)
            {}

            bool operator < (const BoundaryFaceId &other) const
            {
              return adjacent_cell < other.adjacent_cell;
            }
          };

          struct Parameters
          {
            double        CFL_number;
            unsigned int  minimum_substep_number;
            double        diffusion_constant;

            double        stabilization_c_R;
            double        stabilization_alpha;
            double        stabilization_beta;

            bool          smooth_surface_mesh;
            unsigned int  gradient_descent_iterations;
            double        gradient_descent_step_size;
            double        slope_limit;
          };

          Parameters parameters;

          MPI_Comm surface_mpi_comm;

          Triangulation<surface_dim> surface_mesh;

          Triangulation<surface_dim> interface_mesh;

          double global_Omega_diameter;

          std::map<typename Triangulation<surface_dim>::active_cell_iterator, BoundaryFaceId>
          interface_to_volume_mapping;

          std::map<BoundaryFaceId, typename Triangulation<surface_dim>::active_cell_iterator>
          volume_to_interface_mapping;

          std::vector<Point<dim>> interface_coordinates;

          std::unique_ptr<FESystem<surface_dim>> interface_vel_fe;
          DoFHandler<surface_dim>                interface_vel_dof_handler;
          Vector<double>                         interface_velocities;

          std::unique_ptr<FESystem<surface_dim>> surface_vel_fe;
          DoFHandler<surface_dim>                surface_vel_dof_handler;
          Vector<double>                         surface_velocities;

          FESystem<surface_dim>                  coord_fe;
          DoFHandler<surface_dim>                coord_dof_handler;
          Vector<double>                         surface_coordinates;

          FE_Q<surface_dim>                      topo_spatial_fe;
          DoFHandler<surface_dim>                topo_spatial_dof_handler;
          Vector<double>                         topography_spatial;
          Vector<double>                         old_topography_spatial;
          Vector<double>                         old_old_topography_spatial;

          FE_Q<surface_dim>                      topo_material_fe;
          DoFHandler<surface_dim>                topo_material_dof_handler;
          Vector<double>                         topography_material;

          AffineConstraints<double>              convection_diffusion_constraints;
          SparsityPattern                        convection_diffusion_sp;
          SparseMatrix<double>                   convection_diffusion_matrix;
          Vector<double>                         convection_diffusion_rhs;

          AffineConstraints<double>              laplace_beltrami_constraints;
          SparsityPattern                        laplace_beltrami_sp;
          SparseMatrix<double>                   laplace_beltrami_matrix;
          Vector<double>                         laplace_beltrami_rhs;
      };
    }
  }
}

#endif
