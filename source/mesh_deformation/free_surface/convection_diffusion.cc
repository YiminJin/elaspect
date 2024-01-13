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


#include <elaspect/mesh_deformation/free_surface/convection_diffusion.h>
#include <elaspect/mesh_deformation/handler.h>
#include <elaspect/initial_topography/interface.h>
#include <elaspect/boundary_velocity/interface.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/derivative_approximation.h>
#include <deal.II/numerics/data_out.h>

namespace elaspect
{
  namespace MeshDeformation
  {
    namespace FreeSurface
    {
      namespace internal
      {
        namespace GeometryHelper
        {
          static const int child_line_in_quad[4][4]
            = { {  0, -1,  1, -1 }, { -1,  0, -1,  1 },
                {  0,  1, -1, -1 }, { -1, -1,  0,  1 } };

          static const int child_quad_in_hex[6][8]
            = { {  0, -1,  1, -1,  2, -1,  3, -1 },
                { -1,  0, -1,  1, -1,  2, -1,  3 },
                {  0,  2, -1, -1,  1,  5, -1, -1 },
                { -1, -1,  0,  1, -1, -1,  2,  3 },
                {  0,  1,  2,  3, -1, -1, -1, -1 },
                { -1, -1, -1, -1,  0,  1,  2,  3 } };

          static const int quad_vertices_counterclockwise[4][2]
            = { { 2, 0 }, { 1, 3 }, { 0, 1 }, { 3, 2 } };
        }

        /**
         * If a child cell is adjacent to the n-th face of its parent,
         * then its n-th face is also a child of the parent's n-th face.
         * This auxiliary function helps to determine which child of the
         * parent's face the subface is. If the given subface is not
         * at surface of the parent cell, the function will throw an
         * exception in DEBUG mode.
         */
        template <int dim>
        unsigned int child_face_in_cell (const unsigned int face_no,
                                         const unsigned int subcell);

        template <>
        unsigned int child_face_in_cell<2> (const unsigned int face_no,
                                            const unsigned int subcell)
        {
          AssertIndexRange (face_no, 4);
          AssertIndexRange (subcell, 4);
          const int idx = GeometryHelper::child_line_in_quad[face_no][subcell];
          Assert (idx != -1, ExcMessage ("Subcell " 
                                         + Utilities::int_to_string(subcell)
                                         + " is not adjacent to face " 
                                         + Utilities::int_to_string(face_no)
                                         + " in 2d."));

          return idx;
        }


        template <>
        unsigned int child_face_in_cell<3> (const unsigned int face_no,
                                            const unsigned int subcell)
        {
          AssertIndexRange (face_no, 6);
          AssertIndexRange (subcell, 8);
          const int idx = GeometryHelper::child_quad_in_hex[face_no][subcell];
          Assert (idx != -1, ExcMessage ("Subcell " 
                                         + Utilities::int_to_string(subcell)
                                         + " is not adjacent to face " 
                                         + Utilities::int_to_string(face_no)
                                         + " in 3d."));

          return idx;
        }


        /**
         * This function implements a surface subgrid extraction, which is 
         * similar to GridGenerator::extract_boundary_mesh. The differences
         * are:
         * (a) this function only extract the coarsest level of surface mesh;
         * (b) the spacedim of the extracted surface mesh equals dim-1;
         * (c) the boundary ids of the surface mesh consist with the 
         *     volume mesh.
         */
        template <int dim>
        std::map<typename Triangulation<dim>::cell_iterator,
                 typename Triangulation<dim-1>::cell_iterator>
        extract_coarse_surface_mesh (const Triangulation<dim>           &volume_mesh,
                                     Triangulation<dim-1>               &surface_mesh,
                                     const std::set<types::boundary_id> &boundary_ids)
        {
          const unsigned int surface_dim = dim-1;
          surface_mesh.clear();

          // Temporary map from interface cells to volume cells.
          std::vector<typename Triangulation<dim>::cell_iterator> cell_mapping;

          std::vector<Point<surface_dim>>    vertices;
          std::vector<CellData<surface_dim>> cells;
          SubCellData                        subcell_data;
          std::vector<bool> vertex_visited(volume_mesh.n_vertices(), false);

          // Mapping from volume vertex indices to surface vertex indices.
          std::map<unsigned int, unsigned int> vertex_mapping;
          for (const auto &cell : volume_mesh.cell_iterators_on_level(0))
            for (const unsigned int f : cell->face_indices())
            {
              const typename Triangulation<dim>::face_iterator face = cell->face(f);
              if (face->at_boundary() &&
                  boundary_ids.find(face->boundary_id()) != boundary_ids.end())
              {
                CellData<surface_dim> c_data;
                for (const unsigned int v : face->vertex_indices())
                {
                  const unsigned int vertex_index = face->vertex_index(v);
                  if (!vertex_visited[vertex_index])
                  {
                    vertex_visited[vertex_index] = true;
                    const Point<dim> volume_vertex = face->vertex(v);
                    Point<surface_dim> surface_vertex;
                    for (unsigned int  d = 0; d < surface_dim; ++d)
                      surface_vertex[d] = volume_vertex[d];
                    vertices.push_back(surface_vertex);
                    vertex_mapping[vertex_index] = vertices.size() - 1;
                  }

                  c_data.vertices[v] = vertex_mapping[vertex_index];
                  c_data.material_id = cell->material_id();
                  // Discard manifold id.
                }

                // If we start from a 3d mesh, then we need to determin
                // the edges of the surface cell.
                if (dim == 3)
                {
                  // Determine the boundary id of each line bounding the
                  // face. The boundary id of a bounding line is the
                  // boundary id of the other face adjacent to it.
                  std::vector<unsigned int> line_boundary_ids(GeometryInfo<dim>::lines_per_cell,
                                                              numbers::internal_face_boundary_id);
                  for (unsigned int ff = 0; ff < GeometryInfo<dim>::faces_per_cell; ++ff)
                    if (ff != f && ff != GeometryInfo<dim>::opposite_face[f] &&
                        cell->at_boundary(ff))
                    {
                      const types::boundary_id bid = cell->face(ff)->boundary_id();
                      AssertThrow (boundary_ids.find(bid) == boundary_ids.end(),
                                   ExcMessage ("We cannot handle the case that there are cells with "
                                               "more than one face belonging to the free surface."));
                      for (unsigned int l = 0; l < 4; ++l)
                      {
                        const unsigned int line_index = 
                          GeometryInfo<dim>::face_to_cell_lines(ff, l,
                                                                cell->face_orientation(ff),
                                                                cell->face_flip(ff),
                                                                cell->face_rotation(ff));
                        line_boundary_ids[line_index] = cell->face(ff)->boundary_id();
                      }
                    }

                  for (unsigned int l = 0; l < 4; ++l)
                  {
                    // See if we already saw this line from a
                    // neighboring face, either in this or the reverse
                    // orientation. If so, skip it.
                    bool line_found = false;
                    for (unsigned int i=0; i<subcell_data.boundary_lines.size(); ++i)
                      if (((subcell_data.boundary_lines[i].vertices[0]
                            == vertex_mapping[face->line(l)->vertex_index(0)])
                           &&
                           (subcell_data.boundary_lines[i].vertices[1]
                            == vertex_mapping[face->line(l)->vertex_index(1)]))
                          ||
                          ((subcell_data.boundary_lines[i].vertices[0]
                            == vertex_mapping[face->line(l)->vertex_index(1)])
                           &&
                           (subcell_data.boundary_lines[i].vertices[1]
                            == vertex_mapping[face->line(l)->vertex_index(0)])))
                      {
                        line_found = true;
                        break;
                      }

                    if (line_found)
                      continue;

                    CellData<1> line;
                    line.vertices[0] = vertex_mapping[face->line(l)->vertex_index(0)];
                    line.vertices[1] = vertex_mapping[face->line(l)->vertex_index(1)];

                    const unsigned int line_index =
                      GeometryInfo<dim>::face_to_cell_lines (f, l,
                                                             cell->face_orientation(f),
                                                             cell->face_flip(f),
                                                             cell->face_rotation(f));
                    line.boundary_id = line_boundary_ids[line_index];
                    // Discard manifold id.

                    subcell_data.boundary_lines.push_back (line);
                  }
                }

                cells.push_back (c_data);
                cell_mapping.push_back (cell);
              }
            }

          // Create the surface triangulation.
          surface_mesh.create_triangulation (vertices, cells, subcell_data);

          // Return the mapping from volume cells to surface cells.
          std::map<typename Triangulation<dim>::cell_iterator,
                   typename Triangulation<dim-1>::cell_iterator> result_mapping;
          for (const auto &cell : surface_mesh.cell_iterators_on_level(0))
            result_mapping[cell_mapping.at(cell->index())] = cell;

          return result_mapping;
        }


        /**
         * Returns the determinant of a triangle.
         */
        double triangle_determinant (const Point<2> &p1,
                                     const Point<2> &p2,
                                     const Point<2> &p3)
        {
          return (p1[0] - p2[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p2[1]);
        }

    
        DeclException1 (ExcPointOutOfDomain1,
                        Point<1>,
                        << "The program tries to find the active cell around point "
                        << arg1 << ", but the point is not in the surface domain.");

        DeclException1 (ExcPointOutOfDomain2,
                        Point<2>,
                        << "The program tries to find the active cell around point "
                        << arg1 << ", but the point is not in the surface domain.");

        /**
         * Find the active cell that surrounds the given point. Different from the
         * built-in function with the same name in namespace GridTools, this function
         * is based on a recursive strategy, which is less secure but faster, especially
         * when the given point is far away (but not too far) from the hint cell.
         */
        DoFHandler<1>::active_cell_iterator
        find_active_cell_around_point(const Mapping<1>                          &mapping,
                                      const DoFHandler<1>::active_cell_iterator &cell_hint,
                                      const Point<1>                            &point)
        {
          const double tol = 1e-8;
          const Point<1> unit_point = mapping.transform_real_to_unit_cell(cell_hint, point);
          if (unit_point[0] < -tol)
          {
            AssertThrow(!cell_hint->at_boundary(0), ExcPointOutOfDomain1(point));
            DoFHandler<1>::cell_iterator neighbor = cell_hint->neighbor(0);
            if (neighbor->has_children())
              return find_active_cell_around_point(mapping, neighbor->child(1), point);
            else
              return find_active_cell_around_point(mapping, neighbor, point);
          }
          if (unit_point[0] > 1.0 + tol)
          {
            AssertThrow(!cell_hint->at_boundary(1), ExcPointOutOfDomain1(point));
            DoFHandler<1>::cell_iterator neighbor = cell_hint->neighbor(1);
            if (neighbor->has_children())
              return find_active_cell_around_point(mapping, neighbor->child(0), point);
            else
              return find_active_cell_around_point(mapping, neighbor, point);
          }

          return cell_hint;
        }


        DoFHandler<2>::active_cell_iterator
        find_active_cell_around_point(const Mapping<2>                          &mapping,
                                      const DoFHandler<2>::active_cell_iterator &cell_hint,
                                      const Point<2>                            &point)
        {
          const double tol = 1e-8;
          const RefinementCase<2> ref_case (RefinementPossibilities<2>::isotropic_refinement);

          // try transform_real_to_unit_cell first. if the point is inside the present cell
          // then calculations of triangle determinants are unnecessary.
          bool point_inside_cell = false;
          try
          {
            Point<2> unit_point = mapping.transform_real_to_unit_cell(cell_hint, point);
            if (GeometryInfo<2>::is_inside_unit_cell(unit_point, tol))
              point_inside_cell = true;
          }
          catch (...)
          {
            // transforming failed. the point is not inside the present cell
          }

          if (point_inside_cell)
            return cell_hint;

          // determine the position of the point with respect to the present cell
          Point<2> vertices[4];
          for (unsigned int v = 0; v < 4; ++v)
            vertices[v] = mapping.transform_unit_to_real_cell(
              cell_hint, GeometryInfo<2>::unit_cell_vertex(v));

          for (unsigned int i = 0; i < 4; ++i)
          {
            const unsigned int v1 = GeometryHelper::quad_vertices_counterclockwise[i][0],
                               v2 = GeometryHelper::quad_vertices_counterclockwise[i][1];
            if (triangle_determinant(vertices[v1], vertices[v2], point) < tol * (vertices[v1] - vertices[v2]).norm())
            {
              AssertThrow (!cell_hint->at_boundary(i), ExcPointOutOfDomain2(point));
              DoFHandler<2>::cell_iterator neighbor = cell_hint->neighbor(i);
              if (neighbor->has_children())
              {
                unsigned int f = cell_hint->neighbor_of_neighbor(i);
                // simply pick the child adjacent to the 0th subface
                unsigned int c = GeometryInfo<2>::child_cell_on_face(ref_case, f, 0);
                return find_active_cell_around_point(mapping, neighbor->child(c), point);
              }
              else
                return find_active_cell_around_point(mapping, neighbor, point);
            }
          }

          // this should not happen. return something to prevent warning.
          return cell_hint->get_dof_handler().end();
        }


        void
        make_no_normal_flux_constraints(const DoFHandler<1>                   &dof_handler,
                                        const std::vector<types::boundary_id> &boundary_ids,
                                        AffineConstraints<double>             &constraints)
        {
          constraints.clear();
          for (const auto bid : boundary_ids)
            VectorTools::interpolate_boundary_values(dof_handler, 
                                                     bid, 
                                                     Functions::ZeroFunction<1>(1),
                                                     constraints);
          constraints.close();
        }


        void
        make_no_normal_flux_constraints(const DoFHandler<2>                   &dof_handler,
                                        const std::vector<types::boundary_id> &boundary_ids,
                                        AffineConstraints<double>             &constraints)
        {
          std::set<types::boundary_id> bids(boundary_ids.begin(), boundary_ids.end());

          constraints.clear();
          VectorTools::compute_no_normal_flux_constraints (dof_handler,
                                                           /*first vector component = */0,
                                                           bids,
                                                           constraints);
          constraints.close();
        }
      }


      template <int dim>
      ConvectionDiffusion<dim>::ConvectionDiffusion()
        :
        interface_vel_dof_handler(interface_mesh),
        surface_vel_dof_handler(surface_mesh),

        coord_fe(FE_Q<surface_dim>(1), surface_dim),
        coord_dof_handler(surface_mesh),

        topo_spatial_fe(1),
        topo_spatial_dof_handler(surface_mesh),

        topo_material_fe(1),
        topo_material_dof_handler(surface_mesh)
      {}


      template <int dim>
      void ConvectionDiffusion<dim>::create_interface_mesh()
      {
        // Create the coarse mesh.
        const parallel::distributed::Triangulation<dim> &volume_mesh = this->get_triangulation();
        const std::set<types::boundary_id> fs_boundary_ids =
          this->get_mesh_deformation_handler().get_free_surface_boundary_indicators();
        const std::map<typename Triangulation<dim>::cell_iterator, 
                       typename Triangulation<surface_dim>::cell_iterator>
        coarse_cell_mapping = internal::extract_coarse_surface_mesh(volume_mesh,
                                                                    interface_mesh,
                                                                    fs_boundary_ids);

        // Since the free surface might be owned by several mpi processes, we
        // have to gather the information of active cells to the root process.
        const unsigned int cellid_size = sizeof(CellId::binary_type);

        // The number of active cells to be sent should not exceed the number
        // of locally owned active cells.
        std::vector<char> send_data(volume_mesh.n_locally_owned_active_cells() * cellid_size);
        char *send_data_it = &send_data.front();

        for (const auto &cell : volume_mesh.active_cell_iterators())
          if (cell->is_locally_owned() && cell->at_boundary())
            for (const unsigned int face_no : cell->face_indices())
            {
              typename Triangulation<dim>::face_iterator face = cell->face(face_no);
              if (face->at_boundary() &&
                  fs_boundary_ids.find(face->boundary_id()) != fs_boundary_ids.end())
              {
                const CellId::binary_type cellid = cell->id().template to_binary<dim>();
                memcpy(send_data_it, &cellid, cellid_size);
                send_data_it += cellid_size;
              }
            }

        // Gather the data length to all processes.
        const unsigned int global_comm_size = Utilities::MPI::n_mpi_processes (this->get_mpi_communicator());
        int send_data_length = send_data_it - &send_data.front();
        std::vector<int> gathered_data_length (global_comm_size);

        int ierr = MPI_Allgather(&send_data_length, 1, MPI_INT,
                                 &gathered_data_length[0], 1, MPI_INT,
                                 this->get_mpi_communicator());
        AssertThrowMPI(ierr);

        // From the gathered message we can know which processes hold part of 
        // the free surface. For future convenience, create an mpi communicator
        // that contains those processes. 
        std::vector<int> fs_owner_ids;
        for (unsigned int i = 0; i < gathered_data_length.size(); ++i)
          if (gathered_data_length[i] > 0)
            fs_owner_ids.push_back(i);

        // Make sure that the sub-communicator contains process 0.
        if (fs_owner_ids[0] != 0)
          fs_owner_ids.insert(fs_owner_ids.begin(), 0);

        MPI_Group global_mpi_group;
        MPI_Comm_group(this->get_mpi_communicator(), &global_mpi_group);
        MPI_Group surface_mpi_group;
        MPI_Group_incl(global_mpi_group, fs_owner_ids.size(), &fs_owner_ids.front(), &surface_mpi_group);
        MPI_Comm_create(this->get_mpi_communicator(), surface_mpi_group, &surface_mpi_comm);
        
        volume_to_interface_mapping.clear();
        interface_to_volume_mapping.clear();
        // Processes not in surface_mpi_comm can rest now.
        if (surface_mpi_comm == MPI_COMM_NULL)
          return;

        // Calculate the offset of incoming data from each process.
        std::vector<int> offsets(fs_owner_ids.size());
        offsets[0] = 0;
        for (unsigned int i = 1; i < offsets.size(); ++i)
          offsets[i] = offsets[i-1] + gathered_data_length[fs_owner_ids[i-1]];

        std::vector<int> recv_data_length(fs_owner_ids.size());
        for (unsigned int i = 0; i < fs_owner_ids.size(); ++i)
          recv_data_length[i] = gathered_data_length[fs_owner_ids[i]];

        std::vector<char> recv_data(offsets.back() + gathered_data_length[fs_owner_ids.back()]);

        // Gather data to process 0.
        ierr = MPI_Gatherv(&send_data[0], send_data_length, MPI_CHAR,
                           &recv_data[0], &recv_data_length[0], &offsets[0], MPI_CHAR,
                           0, surface_mpi_comm);
        AssertThrowMPI(ierr);

        // Broadcast the information of active cells at surface to all processes
        // in surface_mpi_comm.
        ierr = MPI_Bcast(&recv_data[0], recv_data.size(), MPI_CHAR, 0, surface_mpi_comm);
        AssertThrowMPI(ierr);

        const char *recv_data_it = &recv_data.front();
        const char *const recv_data_end = &recv_data.back() + 1;

        // Record the CellId and MPI rank for each cell for later use.
        std::vector<std::pair<CellId, unsigned int>> cellids;
        unsigned int surface_mpi_rank = 0;
        while (recv_data_it != recv_data_end)
        {
          if (surface_mpi_rank < offsets.size() - 1 &&
              recv_data_it - &recv_data.front() >= offsets[surface_mpi_rank + 1])
            ++surface_mpi_rank;

          CellId::binary_type binary_cellid;
          memcpy (&binary_cellid, recv_data_it, cellid_size);
          cellids.push_back(std::make_pair(CellId(binary_cellid), surface_mpi_rank));
          recv_data_it += cellid_size;
        }

        AssertDimension(surface_mpi_rank, fs_owner_ids.size() - 1);

        // Find out which face is at free surface for each volume cell
        // in the coarse cell map.
        std::map<unsigned int, unsigned int> fs_face_numbers;
        for (auto p = coarse_cell_mapping.begin(); p != coarse_cell_mapping.end(); ++p)
        {
          const typename Triangulation<dim>::cell_iterator cell = p->first;
          for (const unsigned int face_no : cell->face_indices())
            if (cell->at_boundary(face_no))
            {
              const types::boundary_id boundary_id = cell->face(face_no)->boundary_id();
              if (fs_boundary_ids.find(boundary_id) != fs_boundary_ids.end())
              {
                fs_face_numbers[cell->id().get_coarse_cell_id()] = face_no;
                break;
              }
            }
        }

        // Refine the interface mesh.
        const unsigned int max_level = volume_mesh.n_global_levels();
        for (unsigned int level = 0; level < max_level; ++level)
        {
          for (std::vector<std::pair<CellId, unsigned int> >::const_iterator
               it = cellids.begin(); it != cellids.end(); ++it)
          {
            // Get the coarse volume cell.
            const CellId &cellid = it->first;
            const types::coarse_cell_id coarse_cell_id = cellid.get_coarse_cell_id();
            const typename Triangulation<dim>::cell_iterator
            volume_cell = volume_mesh.create_cell_iterator(CellId(coarse_cell_id,
                                                                  std::vector<std::uint8_t>()));
            const unsigned int face_no = fs_face_numbers[coarse_cell_id];

            // Get the coarse interface cell and refine it to the proper level.
            typename Triangulation<surface_dim>::cell_iterator
            interface_cell = coarse_cell_mapping.find(volume_cell)->second;
            const ArrayView<const std::uint8_t> child_indices = cellid.get_child_indices();
            for (unsigned int i = 0; i < child_indices.size(); ++i)
            {
              if (!interface_cell->has_children())
              {
                interface_cell->set_refine_flag();
                break;
              }

              interface_cell = interface_cell->child(
                internal::child_face_in_cell<dim>(face_no, child_indices[i]));
            }
          }

          interface_mesh.execute_coarsening_and_refinement();
        }

        // Build the mapping between volume faces and interface cells.
        for (auto p = cellids.begin(); p != cellids.end(); ++p)
        {
          const CellId &cellid = p->first;
          const types::coarse_cell_id coarse_cell_id = cellid.get_coarse_cell_id();
          const typename Triangulation<dim>::cell_iterator
          volume_cell = volume_mesh.create_cell_iterator(CellId(coarse_cell_id,
                                                                std::vector<std::uint8_t>()));
          const unsigned int face_no = fs_face_numbers[coarse_cell_id];

          typename Triangulation<surface_dim>::cell_iterator
          interface_cell = coarse_cell_mapping.find(volume_cell)->second;
          const ArrayView<const std::uint8_t> child_indices = cellid.get_child_indices();
          for (unsigned int i = 0; i < child_indices.size(); ++i)
            interface_cell = interface_cell->child(
              internal::child_face_in_cell<dim>(face_no, child_indices[i]));

          BoundaryFaceId faceid(cellid, face_no, p->second);
          interface_to_volume_mapping[interface_cell] = faceid;
          // volume_to_interface_mapping only stores volume faces that are
          // locally owned.
          if (p->second == Utilities::MPI::this_mpi_process(surface_mpi_comm))
            volume_to_interface_mapping[faceid] = interface_cell;
        }

        AssertDimension(interface_to_volume_mapping.size(), interface_mesh.n_active_cells());

        // Adjust the coordinates of the interface mesh.
        std::vector<bool> vertex_visited(interface_mesh.n_vertices(), false);
        for (const auto &interface_cell : interface_mesh.active_cell_iterators())
        {
          // Get the corresponding cell in surface mesh (which is not necessarily active).
          typename Triangulation<surface_dim>::cell_iterator
          surface_cell = surface_mesh.create_cell_iterator(interface_cell->id());

          for (const unsigned int v : interface_cell->vertex_indices())
            if (!vertex_visited[interface_cell->vertex_index(v)])
            {
              vertex_visited[interface_cell->vertex_index(v)] = true;

              // Find the active cell in surface mesh that contains the target vertex.
              typename DoFHandler<surface_dim>::cell_iterator dof_cell(&surface_mesh, 
                                                                       surface_cell->level(), 
                                                                       surface_cell->index(), 
                                                                       &coord_dof_handler);
              while (dof_cell->has_children())
                dof_cell = dof_cell->child(v);

              Point<surface_dim> point;
              for (unsigned int d = 0; d < surface_dim; ++d)
                point[d] = surface_coordinates[dof_cell->vertex_dof_index(v, d)];
              interface_cell->vertex(v) = point;
            }
        }

        // Finally, re-initialize the containers associated with the interface mesh.
        interface_vel_dof_handler.distribute_dofs(*interface_vel_fe);
        interface_velocities.reinit(interface_vel_dof_handler.n_dofs());
        interface_coordinates.resize(interface_mesh.n_vertices());
      }


      template <int dim>
      void
      ConvectionDiffusion<dim>::gather_data_to_root_process()
      {
        // This function is executed only on processes that hold
        // part of the free surface.
        if (surface_mpi_comm == MPI_COMM_NULL)
          return;

        const Triangulation<dim> &volume_mesh              = this->get_triangulation();
        const FiniteElement<dim> &volume_fe                = this->get_fe();
        const DoFHandler<dim> &volume_dof_handler          = this->get_dof_handler();
        const TrilinosWrappers::MPI::BlockVector &solution = this->get_solution();

        const double time_step = this->get_timestep();

        // Only the velocities (incremental displacements divided by time step) 
        // ought to be passed to interface, so we need to know which dof indices 
        // belong to the displacement components.
        const FiniteElement<dim> &disp_fe = 
          volume_fe.base_element(this->introspection().base_elements.displacement);
        const unsigned int dofs_per_vertex = disp_fe.dofs_per_vertex * dim,
                           dofs_per_line   = disp_fe.dofs_per_line * dim,
                           dofs_per_quad   = disp_fe.dofs_per_quad * dim;

        // Arrays storing the global dof indices of vertices, lines and quads.
        std::vector<types::global_dof_index> vertex_dof_indices(dofs_per_vertex),
                                             line_dof_indices(dofs_per_line),
                                             quad_dof_indices(dofs_per_quad);

        // Arrays storing the dof values of vertices, lines and quads.
        std::vector<double> vertex_dof_values(dofs_per_vertex),
                            line_dof_values(dofs_per_line),
                            quad_dof_values(dofs_per_quad);

        // Arrays storing the dof indices and values of the data to be
        // transported in one cell.
        std::vector<types::global_dof_index> indices;
        std::vector<double> values;

        // Flag all the vertices (and lines in 3D) that has been visited to
        // make sure that each dof on the surface is recorded only once.
        std::vector<bool> vertex_user_flags(volume_mesh.n_vertices(), false);
        if (dim == 3)
        {
          for (auto p = volume_to_interface_mapping.begin(); p != volume_to_interface_mapping.end(); ++p)
          {
            const typename Triangulation<dim>::cell_iterator 
            cell = volume_mesh.create_cell_iterator(p->first.adjacent_cell);
            Assert (cell->is_locally_owned() && cell->is_active(), ExcInternalError());

            typename Triangulation<dim>::face_iterator face = cell->face(p->first.face_no);
            for (unsigned int l = 0; l < 4; ++l)
              face->line(l)->clear_user_flag();
          }
        }

        // Data arrays to be sent and received.
        std::vector<char> send_data, recv_data;

        // Since not all dof values are sent for each cell, we 
        // need a signal to remind the receiver which dofs have 
        // been sent. The signal should be able to indicate 
        // whether the dof values on each vertex, line and quad
        // (in 3D) have been sent, so a digital type composed of
        // 16 bits is sufficient.
        const unsigned int signal_size = sizeof(std::uint16_t);

        // All the data ought to be sent for each cell includes
        // a signal, two integers indicating the level and index,
        // and some float numbers indicating the velocity dof values.
        // The total number of float numbers equals to the number of
        // dofs of the interface dof handler.
        const unsigned int total_length
          = interface_mesh.n_active_cells() * (signal_size + sizeof(unsigned int) * 2) +
            interface_vel_dof_handler.n_dofs() * sizeof(double);

        send_data.resize(total_length);
        void *send_data_it = static_cast<void *>(&send_data.front());
        for (auto p = volume_to_interface_mapping.begin(); p != volume_to_interface_mapping.end(); ++p)
        {
          typename DoFHandler<dim>::cell_iterator 
          cell(*(volume_mesh.create_cell_iterator(p->first.adjacent_cell)), &volume_dof_handler);
          typename DoFHandler<dim>::face_iterator face = cell->face(p->first.face_no);

          indices.clear();
          std::uint16_t signal = 0;
          unsigned int n_bits = 0;
          // Get dof indices on vertices if they have not been recorded yet.
          for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face; ++v, ++n_bits)
            if (!vertex_user_flags[face->vertex_index(v)])
            {
              vertex_user_flags[face->vertex_index(v)] = true;

              signal |= (1<<n_bits);
              for (unsigned int i = 0; i < dofs_per_vertex; ++i)
                vertex_dof_indices[i] = face->vertex_dof_index(v, i);
              indices.insert(indices.end(), vertex_dof_indices.begin(), vertex_dof_indices.end());
            }

          // If the interface mesh is in 1D, then each cell is a line. Get 
          // the dof indices of the line.
          if (dim == 2)
          {
            signal |= (1<<n_bits);
            for (unsigned int i = 0; i < dofs_per_line; ++i)
              line_dof_indices[i] = face->dof_index(i);
            indices.insert(indices.end(), line_dof_indices.begin(), line_dof_indices.end());
          }
          // If the interface mesh is in 2D, then each cell is a quadrilateral.
          // Get the dof indices of the quadrilateral and the four lines
          // bounding it, if they have not been recorded yet.
          else
          {
            for (unsigned int l = 0; l < 4; ++l, ++n_bits)
            {
              auto line = face->line(l);
              if (!line->user_flag_set())
              {
                line->set_user_flag();
                signal |= (1<<n_bits);
                for (unsigned int i = 0; i < dofs_per_line; ++i)
                  line_dof_indices[i] = line->dof_index(i);
                indices.insert(indices.end(), line_dof_indices.begin(), line_dof_indices.end());
              }
            }

            signal |= (1<<n_bits);
            for (unsigned int i=0; i<dofs_per_quad; ++i)
              quad_dof_indices[i] = face->dof_index(i);
            indices.insert(indices.end(), quad_dof_indices.begin(), quad_dof_indices.end());
          }

          // Write the signal.
          memcpy(send_data_it, &signal, signal_size);
          send_data_it = static_cast<char *>(send_data_it) + signal_size;
          // Write the level and index of the interface cell.
          unsigned int *cell_data = static_cast<unsigned int *>(send_data_it);
          *cell_data++ = p->second->level();
          *cell_data++ = p->second->index();
          // Write the displacement increments.
          values.resize(indices.size());
          solution.extract_subvector_to(indices, values);
          double *disp_data = reinterpret_cast<double *>(cell_data);
          for (unsigned int i = 0; i < values.size(); ++i)
            *disp_data++ = values[i];
          send_data_it = static_cast<void *>(disp_data);
        }

        // Gather data length to the root process.
        const unsigned int fs_comm_size = Utilities::MPI::n_mpi_processes(surface_mpi_comm);
        unsigned int send_data_length = static_cast<char *>(send_data_it) - &send_data[0];
        std::vector<int> gathered_data_length(fs_comm_size);
        int ierr = MPI_Gather (&send_data_length, 1, MPI_INT,
                               &gathered_data_length[0], 1, MPI_INT,
                               0, surface_mpi_comm);
        AssertThrowMPI(ierr);

        // Calculate the offset of incoming data from each process.
        std::vector<int> offsets(fs_comm_size);
        offsets[0] = 0;
        for (unsigned int i = 1; i < fs_comm_size; ++i)
          offsets[i] = offsets[i-1] + gathered_data_length[i-1];

        // Note that the total length of data to be received does not equal
        // to total_length we calculated above, because vertices and lines
        // shared by different subdomains are visited more than once.
        recv_data.resize (offsets.back() + gathered_data_length.back());

        // Gather data to the root process.
        ierr = MPI_Gatherv (&send_data[0], send_data_length, MPI_CHAR,
                            &recv_data[0], &gathered_data_length[0], &offsets[0], MPI_CHAR,
                            0, surface_mpi_comm);
        AssertThrowMPI(ierr);

        // All the processes except for the root process can rest now.
        if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) != 0)
          return;

        // Read the received data.
        const void *recv_data_it = static_cast<const void *>(&recv_data.front());
        const void *const recv_data_end = static_cast<const void *>(&recv_data.back() + 1);

        std::uint16_t signal;
        while (recv_data_it != recv_data_end)
        {
          // Read the signal.
          memcpy (&signal, recv_data_it, signal_size);
          recv_data_it = static_cast<const char *>(recv_data_it) + signal_size;
          // Read the level and index of the interface cell.
          const unsigned int *cell_data = static_cast<const unsigned int *>(recv_data_it);
          const unsigned int level = *cell_data++;
          const unsigned int index = *cell_data++;
          const typename DoFHandler<surface_dim>::active_cell_iterator
          cell(&interface_mesh, level, index, &interface_vel_dof_handler);
          // Read the displacement increments.
          const double *disp_data = reinterpret_cast<const double *>(cell_data);
          unsigned int n_bits = 0;
          indices.clear();
          values.clear();
          for (unsigned int v = 0; v < GeometryInfo<surface_dim>::vertices_per_cell; ++v, ++n_bits)
            if (signal & (1<<n_bits))
            {
              for (unsigned int i = 0; i < dofs_per_vertex; ++i)
              {
                vertex_dof_indices[i] = cell->vertex_dof_index(v, i);
                vertex_dof_values[i] = *disp_data++;
              }
              indices.insert(indices.end(), vertex_dof_indices.begin(), vertex_dof_indices.end());
              values.insert(values.end(), vertex_dof_values.begin(), vertex_dof_values.end());
            }
          if (dim == 2)
          {
            Assert(signal & (1<<n_bits), ExcInternalError());
            for (unsigned int i = 0; i < dofs_per_line; ++i)
            {
              line_dof_indices[i] = cell->dof_index(i);
              line_dof_values[i] = *disp_data++;
            }
            indices.insert(indices.end(), line_dof_indices.begin(), line_dof_indices.end());
            values.insert(values.end(), line_dof_values.begin(), line_dof_values.end());
          }
          else
          {
            for (unsigned int l = 0; l < 4; ++l, ++n_bits)
              if (signal & (1<<n_bits))
              {
                auto line = cell->line(l);
                for (unsigned int i = 0; i < dofs_per_line; ++i)
                {
                  line_dof_indices[i] = line->dof_index(i);
                  line_dof_values[i] = *disp_data++;
                }
                indices.insert(indices.end(), line_dof_indices.begin(), line_dof_indices.end());
                values.insert(values.end(), line_dof_values.begin(), line_dof_values.end());
              }

            Assert (signal & (1<<n_bits), ExcInternalError());
            for (unsigned int i = 0; i < dofs_per_quad; ++i)
            {
              quad_dof_indices[i] = cell->dof_index(i);
              quad_dof_values[i] = *disp_data++;
            }
            indices.insert(indices.end(), quad_dof_indices.begin(), quad_dof_indices.end());
            values.insert(values.end(), quad_dof_values.begin(), quad_dof_values.end());
          }

          // Store the velocity values.
          for (unsigned int i = 0; i < indices.size(); ++i)
            interface_velocities[indices[i]] = values[i] / time_step;

          recv_data_it = static_cast<const void *>(disp_data);
        }
      }


      template <int dim>
      void
      ConvectionDiffusion<dim>::project_velocities_onto_surface()
      {
        // the projection is executed on process 0.
        if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) != 0) 
          return;

        // set up the mass matrix.
        DynamicSparsityPattern dsp(surface_vel_dof_handler.n_dofs(),
                                   surface_vel_dof_handler.n_dofs());
        DoFTools::make_sparsity_pattern(surface_vel_dof_handler, dsp);

        SparsityPattern sparsity_pattern;
        sparsity_pattern.copy_from(dsp);

        SparseMatrix<double> mass_matrix(sparsity_pattern);

        // set up the right-hand side vector
        Vector<double> rhs(surface_vel_dof_handler.n_dofs());

        // stuff for assembling
        QGauss<surface_dim> quadrature(surface_vel_fe->degree + 1);
        FEValues<surface_dim> fe_values(*surface_vel_fe, quadrature,
                                        update_values |
                                        update_quadrature_points |
                                        update_JxW_values);

        const unsigned int n_q_points              = quadrature.size();
        const unsigned int surface_dofs_per_cell   = surface_vel_fe->dofs_per_cell;
        const unsigned int interface_dofs_per_cell = interface_vel_fe->dofs_per_cell;

        FullMatrix<double> cell_matrix(surface_dofs_per_cell, surface_dofs_per_cell);
        Vector<double> cell_rhs(surface_dofs_per_cell);

        std::vector<types::global_dof_index> interface_dof_indices(interface_dofs_per_cell);

        for (const auto &surface_cell : surface_vel_dof_handler.active_cell_iterators())
        {
          // take the corresponding interface cell as hint cell. note that the 
          // interface cell may be coarser than the surface cell, so we cannot
          // use Triangulation<dim>::create_cell_iterator() directly.
          CellId cellid = surface_cell->id();
          const ArrayView<const std::uint8_t> child_indices = cellid.get_child_indices();
          const typename Triangulation<surface_dim>::cell_iterator
          interface_cell = interface_mesh.create_cell_iterator(CellId(cellid.get_coarse_cell_id(),
                                                                      std::vector<std::uint8_t>()));
          typename DoFHandler<surface_dim>::cell_iterator
          cell_hint(*interface_cell, &interface_vel_dof_handler);
          unsigned int level = 0;
          while (cell_hint->has_children())
            cell_hint = cell_hint->child(child_indices[level++]);

          fe_values.reinit(surface_cell);

          cell_matrix = 0;
          cell_rhs = 0;
          for (unsigned int q = 0; q < n_q_points; ++q)
          {
            // find the interface cell that covers the quadrature point
            const Point<surface_dim> &q_point = fe_values.quadrature_point(q);
            typename DoFHandler<surface_dim>::active_cell_iterator
            interface_cell = internal::find_active_cell_around_point(
              StaticMappingQ1<surface_dim>::mapping, cell_hint, q_point);
            interface_cell->get_dof_indices(interface_dof_indices);

            const Point<surface_dim> unit_point = StaticMappingQ1<surface_dim>::mapping
              .transform_real_to_unit_cell (interface_cell, q_point);

            // compute the velocity at the quadrature point
            Tensor<1,dim> velocity;
            for (unsigned int i = 0; i < interface_dofs_per_cell; ++i)
            {
              const unsigned int comp = interface_vel_fe->system_to_component_index(i).first;
              velocity[comp] += interface_vel_fe->shape_value(i, unit_point)
                                * interface_velocities[interface_dof_indices[i]];
            }

            for (unsigned int i=0; i<surface_dofs_per_cell; ++i)
            {
              const unsigned int comp_i = surface_vel_fe->system_to_component_index(i).first;
              cell_rhs(i) += fe_values.shape_value(i,q) * velocity[comp_i] * fe_values.JxW(q);

              for (unsigned int j = 0; j < surface_dofs_per_cell; ++j)
              {
                const unsigned int comp_j = surface_vel_fe->system_to_component_index(j).first;
                if (comp_i == comp_j)
                  cell_matrix(i,j) += fe_values.shape_value(i,q) * fe_values.shape_value(j,q) * fe_values.JxW(q);
              }
            }
          }

          surface_cell->distribute_local_to_global(cell_matrix, cell_rhs, mass_matrix, rhs);
        }

        // Solve the linear system.
        SolverControl solver_control(rhs.size(), 1e-8 * rhs.l2_norm());
        SolverCG<Vector<double>> cg(solver_control);

        PreconditionSSOR<SparseMatrix<double>> preconditioner;
        preconditioner.initialize(mass_matrix, 1.2);

        cg.solve(mass_matrix, surface_velocities, rhs, preconditioner);
      }


      template <int dim>
      double 
      ConvectionDiffusion<dim>::
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
                                   const double                              cell_diameter) const
      {
        if (global_u_infty < 1e-50)
          return 0;

        const unsigned int n_q_points = old_z_values.size();

        double max_residual = 0;
        double max_velocity = 0;

        for (unsigned int q = 0; q < n_q_points; ++q)
        {
          const double u_grad_z = velocity_values[q] * old_z_grads[q];

          const double dz_dt = (old_z_values[q] - old_old_z_values[q]) / step_size;

          const double kappa_Delta_z = parameters.diffusion_constant * old_z_laplacian;

          double residual = std::abs(dz_dt + u_grad_z - kappa_Delta_z - uplift_rates[q]);
          if (parameters.stabilization_alpha == 2)
            residual *= std::abs(old_z_values[q] - global_z_average);

          max_residual = std::max(residual, max_residual);
          max_velocity = std::max(velocity_values[q].norm(), max_velocity);
        }

        const double max_viscosity = (parameters.stabilization_beta * max_velocity * cell_diameter);
        if (std::abs(global_entropy_variation) < 1e-50 ||
            std::abs(global_z_variation) < 1e-50)
          return max_viscosity;

        double entropy_viscosity;
        if (parameters.stabilization_alpha == 2)
          entropy_viscosity = (parameters.stabilization_c_R *
                               cell_diameter * cell_diameter *
                               max_residual /
                               global_entropy_variation);
        else
          entropy_viscosity = (parameters.stabilization_c_R *
                               cell_diameter * global_Omega_diameter *
                               max_velocity * max_residual /
                               (global_u_infty * global_z_variation));

        return std::min(max_viscosity, entropy_viscosity);
      }


      template <int dim>
      void 
      ConvectionDiffusion<dim>::compute_topographic_evolution()
      {
        if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) != 0)
          return;

        const unsigned int n_substeps = subdivide_time_step();
        const double substep_size = this->get_timestep() / n_substeps;

        // stuff for assembly
        const QGauss<surface_dim> quadrature(2);
        FEValues<surface_dim> topo_fe_values(topo_spatial_fe, quadrature,
                                             update_values | update_gradients |
                                             update_hessians | update_JxW_values);
        FEValues<surface_dim> vel_fe_values(*surface_vel_fe, quadrature, update_values);

        FEValuesExtractors::Vector u_extractor(0);
        FEValuesExtractors::Scalar r_extractor(surface_dim);

        const unsigned int n_q_points    = topo_fe_values.n_quadrature_points;
        const unsigned int dofs_per_cell = topo_fe_values.dofs_per_cell;

        std::vector<double> old_z_values(n_q_points), old_old_z_values(n_q_points);
        std::vector<Tensor<1,surface_dim>> old_z_grads(n_q_points);
        std::vector<Tensor<1,surface_dim>> velocity_values(n_q_points);
        std::vector<double> uplift_rates(n_q_points);

        std::vector<double> phi(dofs_per_cell);
        std::vector<Tensor<1,surface_dim>> grad_phi(dofs_per_cell);

        std::vector<types::global_dof_index> cell_dof_indices(dofs_per_cell);

        FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
        Vector<double> cell_rhs(dofs_per_cell);

        for (unsigned int substep_number = 0; substep_number < n_substeps; ++substep_number)
        {
          // get the global range of topography
          std::pair<double, double> global_z_range;
          global_z_range.first =  *std::min_element(old_topography_spatial.begin(),
                                                    old_topography_spatial.end());
          global_z_range.second = *std::max_element(old_topography_spatial.begin(),
                                                    old_topography_spatial.end());

          const double global_z_average = (global_z_range.first + global_z_range.second) / 2;

          // get the entropy variation
          double global_entropy_variation = numbers::signaling_nan<double>();
          if (parameters.stabilization_alpha == 2)
          {
            double min_entropy = std::numeric_limits<double>::max(),
                   max_entropy = -std::numeric_limits<double>::max(),
                   area = 0,
                   entropy_integrated = 0;

            for (const auto &topo_cell : topo_spatial_dof_handler.active_cell_iterators())
            {
              topo_fe_values.reinit(topo_cell);
              topo_fe_values.get_function_values(old_topography_spatial, old_z_values);
              for (unsigned int q = 0; q < n_q_points; ++q)
              {
                const double entropy = (old_z_values[q] - global_z_average) *
                                       (old_z_values[q] - global_z_average);

                min_entropy = std::min(min_entropy, entropy);
                max_entropy = std::max(max_entropy, entropy);
                area += topo_fe_values.JxW(q);
                entropy_integrated += topo_fe_values.JxW(q) * entropy;
              }
            }

            const double average_entropy = entropy_integrated / area;
            global_entropy_variation = std::max(max_entropy - average_entropy,
                                                average_entropy - min_entropy);
          }

          // get the maximal velocity
          double global_max_velocity = 0;
          for (const auto &vel_cell : surface_vel_dof_handler.active_cell_iterators())
          {
            vel_fe_values.reinit(vel_cell);
            vel_fe_values[u_extractor].get_function_values(surface_velocities, velocity_values);
            
            for (unsigned int q = 0; q < n_q_points; ++q)
              global_max_velocity = std::max(global_max_velocity, velocity_values[q].norm());
          }

          convection_diffusion_matrix = 0;
          convection_diffusion_rhs = 0;

          typename DoFHandler<surface_dim>::active_cell_iterator
          topo_cell = topo_spatial_dof_handler.begin_active(),
          topo_endc = topo_spatial_dof_handler.end();
          typename DoFHandler<surface_dim>::active_cell_iterator
          vel_cell = surface_vel_dof_handler.begin_active();
          for (; topo_cell != topo_endc; ++topo_cell, ++vel_cell)
          {
            topo_fe_values.reinit(topo_cell);
            topo_fe_values.get_function_values(old_topography_spatial, old_z_values);
            topo_fe_values.get_function_values(old_old_topography_spatial, old_old_z_values);
            topo_fe_values.get_function_gradients(old_topography_spatial, old_z_grads);
            topo_cell->get_dof_indices(cell_dof_indices);
            
            vel_fe_values.reinit(vel_cell);
            vel_fe_values[u_extractor].get_function_values(surface_velocities, velocity_values);
            vel_fe_values[r_extractor].get_function_values(surface_velocities, uplift_rates);

            // compute the approximated laplacian of surface topography. we do not 
            // use function get_function_laplacians() because in piecewise bilinear
            // polynomial space the laplacian of shape function is always 0. note that
            // the laplacian is used for computing the residual, which is necessary
            // only when the diffusion constant is not 0.
            Tensor<2,surface_dim> old_z_hessian;
            if (parameters.diffusion_constant > 0)
              DerivativeApproximation::approximate_derivative_tensor(topo_spatial_dof_handler,
                                                                     old_topography_spatial,
                                                                     topo_cell,
                                                                     old_z_hessian);
            const double old_z_laplacian = trace(old_z_hessian);

            // compute the artificial viscosity
            const double artificial_viscosity =
              compute_artificial_viscosity(old_z_values,
                                           old_old_z_values,
                                           old_z_grads,
                                           old_z_laplacian,
                                           velocity_values,
                                           uplift_rates,
                                           global_max_velocity,
                                           global_z_range.second - global_z_range.first,
                                           global_z_average,
                                           global_entropy_variation,
                                           substep_size,
                                           topo_cell->diameter());

            cell_matrix = 0;
            cell_rhs = 0;
            for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
              {
                phi[k]      = topo_fe_values.shape_value(k, q);
                grad_phi[k] = topo_fe_values.shape_grad(k, q);
              }

              const double diffusion_constant = std::max(parameters.diffusion_constant,
                                                         artificial_viscosity);
              const double JxW = topo_fe_values.JxW(q);

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                cell_rhs(i)
                += phi[i] * (old_z_values[q] + substep_size * uplift_rates[q]) * JxW;

                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  cell_matrix(i, j)
                  += (phi[i] * phi[j]
                      +
                      substep_size * 
                      (phi[i] * (velocity_values[q] * grad_phi[j]) +
                       diffusion_constant * (grad_phi[i] * grad_phi[j]))
                     )
                     * JxW;
              }
            }

            convection_diffusion_constraints.distribute_local_to_global(cell_matrix, cell_rhs,
                                                                        cell_dof_indices,
                                                                        convection_diffusion_matrix,
                                                                        convection_diffusion_rhs);
          }

          // Solve the linear system with GMRES solver.
          const double rhs_norm = convection_diffusion_rhs.l2_norm();
          
          // Check if the right-hand side is zero.
          if (rhs_norm <= std::numeric_limits<double>::min())
          {
            topography_spatial = 0.0;
            continue;
          }

          SolverControl solver_control(1000, 1e-8 * rhs_norm);
          SolverGMRES<Vector<double>> solver(solver_control);

          convection_diffusion_constraints.set_zero(topography_spatial);

          solver.solve(convection_diffusion_matrix,
                       topography_spatial,
                       convection_diffusion_rhs,
                       PreconditionIdentity());

          convection_diffusion_constraints.distribute(topography_spatial);

          old_old_topography_spatial = old_topography_spatial;
          old_topography_spatial = topography_spatial;
        }
      }


      template <int dim>
      unsigned int
      ConvectionDiffusion<dim>::subdivide_time_step() const
      {
        const QIterated<surface_dim> quadrature_formula(QTrapezoid<1>(),
                                                        surface_vel_fe->degree);
        FEValues<surface_dim> fe_values(*surface_vel_fe, quadrature_formula, update_values);
        FEValuesExtractors::Vector u_extractor(0);

        const unsigned int n_q_points = quadrature_formula.size();
        std::vector<Tensor<1,surface_dim>> velocity_values(n_q_points);

        double max_speed_over_meshsize = 0;
        double min_meshsize = std::numeric_limits<double>::max();

        typename DoFHandler<surface_dim>::active_cell_iterator
        vel_cell = surface_vel_dof_handler.begin_active(),
        vel_endc = surface_vel_dof_handler.end();
        typename DoFHandler<surface_dim>::active_cell_iterator
        topo_cell = topo_spatial_dof_handler.begin_active();
        for (; vel_cell != vel_endc; ++vel_cell, ++topo_cell)
        {
          fe_values.reinit(vel_cell);
          fe_values[u_extractor].get_function_values(surface_velocities, velocity_values);

          const double meshsize = vel_cell->minimum_vertex_distance();
          min_meshsize = std::min(min_meshsize, meshsize);

          double u_norm = 0;
          for (unsigned int q = 0; q < n_q_points; ++q)
            u_norm = std::max(velocity_values[q].norm(), u_norm);

          max_speed_over_meshsize = std::max(max_speed_over_meshsize,
                                             u_norm / meshsize);
        }

        double min_convection_timestep = std::numeric_limits<double>::max();
        double min_conduction_timestep = std::numeric_limits<double>::max();

        if (max_speed_over_meshsize > 0)
          min_convection_timestep = parameters.CFL_number
                                    / (surface_vel_fe->degree * max_speed_over_meshsize);

        if (parameters.diffusion_constant > 0)
          min_conduction_timestep = parameters.CFL_number * pow(min_meshsize, 2)
                                    / parameters.diffusion_constant;

        double new_time_step = std::min(min_convection_timestep, 
                                        min_conduction_timestep);

        unsigned int n_subdivisions = std::ceil(this->get_timestep() / new_time_step);

        return std::max(parameters.minimum_substep_number, n_subdivisions);
      }


      template <int dim>
      void
      ConvectionDiffusion<dim>::execute_mesh_smoothing()
      {
        if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) != 0)
          return;

        // Update the vertex coordinates of surface mesh.
        std::vector<bool> vertex_visited(surface_mesh.n_vertices(), false);
        for (const auto &cell : coord_dof_handler.active_cell_iterators())
          for (const unsigned int v : cell->vertex_indices())
            if (!vertex_visited[cell->vertex_index(v)])
            {
              vertex_visited[cell->vertex_index(v)] = true;

              const Point<surface_dim> &point = cell->vertex(v);
              for (unsigned int d = 0; d < surface_dim; ++d)
                surface_coordinates[cell->vertex_dof_index(v, d)] = point[d];
            }

        // stuff for assembly
        QGauss<surface_dim> quadrature(2);
        FEValues<surface_dim> xy_fe_values(coord_fe, quadrature,
                                           update_values | update_gradients | update_JxW_values);
        FEValues<surface_dim> z_fe_values(topo_material_fe, quadrature, update_gradients);

        const FEValuesExtractors::Vector xy_extractor(0);

        const unsigned int n_q_points    = xy_fe_values.n_quadrature_points;
        const unsigned int dofs_per_cell = xy_fe_values.dofs_per_cell;

        std::vector<Tensor<1,surface_dim>> phi(dofs_per_cell);
        std::vector<Tensor<2,surface_dim>> Grad_phi(dofs_per_cell);

        std::vector<Tensor<2,surface_dim>> Grad_xy(n_q_points);
        std::vector<Tensor<1,surface_dim>> Grad_z(n_q_points);

        FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
        Vector<double> cell_rhs(dofs_per_cell);
        std::vector<types::global_dof_index> cell_dof_indices (dofs_per_cell);

        Vector<double> solution(coord_dof_handler.n_dofs());

        // Execute the Gradient Descent algorithm.
        for (unsigned int iter = 0; iter < parameters.gradient_descent_iterations; ++iter)
        {
          update_topography_in_material_coordinates();

          laplace_beltrami_matrix = 0;
          laplace_beltrami_rhs = 0;

          typename DoFHandler<surface_dim>::active_cell_iterator
          xy_cell = coord_dof_handler.begin_active(),
          xy_endc = coord_dof_handler.end();
          typename DoFHandler<surface_dim>::active_cell_iterator
          z_cell = topo_material_dof_handler.begin_active();
          for (; xy_cell != xy_endc; ++xy_cell, ++z_cell)
          {
            xy_fe_values.reinit(xy_cell);
            xy_fe_values[xy_extractor].get_function_gradients(surface_coordinates, Grad_xy);
            xy_cell->get_dof_indices(cell_dof_indices);

            z_fe_values.reinit(z_cell);
            z_fe_values.get_function_gradients(topography_material, Grad_z);

            cell_matrix = 0;
            cell_rhs = 0;
            for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
              {
                phi[k]      = xy_fe_values[xy_extractor].value(k, q);
                Grad_phi[k] = xy_fe_values[xy_extractor].gradient(k, q);
              }

              double slope  = Grad_z[q].norm();
              double factor = 1. / (1. + slope / parameters.slope_limit);
              Tensor<1,surface_dim> Grad_z_star = Grad_z[q] * factor;

              SymmetricTensor<2,surface_dim> G;
              for (unsigned int i = 0; i < surface_dim; ++i)
                for (unsigned int j = i; j < surface_dim; ++j)
                {
                  for (unsigned int d = 0; d < surface_dim; ++d)
                    G[i][j] += Grad_xy[q][d][i] * Grad_xy[q][d][j];

                  G[i][j] += Grad_z_star[i] * Grad_z_star[j];
                }

              double g = 1. / determinant(G);
              G *= std::sqrt(g);

              const double JxW = xy_fe_values.JxW(q);

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                cell_rhs(i) -= (G * symmetrize(transpose(Grad_xy[q]) * Grad_phi[i])) * JxW;

                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  cell_matrix(i,j) += (G * symmetrize(transpose(Grad_phi[j]) * Grad_phi[i])) * JxW;
              }
            }

            laplace_beltrami_constraints.distribute_local_to_global(cell_matrix, cell_rhs,
                                                                    cell_dof_indices,
                                                                    laplace_beltrami_matrix,
                                                                    laplace_beltrami_rhs);             
          }

          // Solve the linear system with CG solver and SSOR preconditioner.
          const double rhs_norm = laplace_beltrami_rhs.l2_norm();

          // Check if the right-hand side is zero.
          if (rhs_norm < std::numeric_limits<double>::min())
          {
            surface_coordinates = 0.0;
            continue;
          }

          PreconditionSSOR<SparseMatrix<double>> preconditioner;
          preconditioner.initialize(laplace_beltrami_matrix, 1.2);

          SolverControl solver_control(laplace_beltrami_rhs.size(), 1e-8 * rhs_norm);
          SolverCG<Vector<double>> solver(solver_control);

          solution = 0;
          solver.solve(laplace_beltrami_matrix,
                       solution,
                       laplace_beltrami_rhs,
                       preconditioner);

          laplace_beltrami_constraints.distribute(solution);

          // Update surface coordinates.
          surface_coordinates.add(parameters.gradient_descent_step_size, solution);
        }

        update_topography_in_material_coordinates();
      }


      template <int dim>
      void
      ConvectionDiffusion<dim>::update_topography_in_material_coordinates()
      {
        const unsigned int dofs_per_cell = topo_spatial_fe.dofs_per_cell;
        Vector<double> topo_dof_values(dofs_per_cell);

        std::vector<bool> vertex_visited(surface_mesh.n_vertices(), false);

        typename DoFHandler<surface_dim>::active_cell_iterator
        xy_cell = coord_dof_handler.begin_active(),
        xy_endc = coord_dof_handler.end();
        typename DoFHandler<surface_dim>::active_cell_iterator
        z_cell = topo_material_dof_handler.begin_active();
        for (; xy_cell != xy_endc; ++xy_cell, ++z_cell)
        {
          for (const unsigned int v : xy_cell->vertex_indices())
            if (!vertex_visited[xy_cell->vertex_index(v)])
            {
              vertex_visited[xy_cell->vertex_index(v)] = true;

              Point<surface_dim> point;
              for (unsigned int d = 0; d < surface_dim; ++d)
                point[d] = surface_coordinates[xy_cell->vertex_dof_index(v, d)];

              // Find the cell that covers the vertex in spatial coordinates.
              typename DoFHandler<surface_dim>::active_cell_iterator
              cell_hint(&surface_mesh, xy_cell->level(), xy_cell->index(), &topo_spatial_dof_handler);

              typename DoFHandler<surface_dim>::active_cell_iterator
              spatial_cell = internal::find_active_cell_around_point(
                StaticMappingQ1<surface_dim>::mapping, cell_hint, point);

              Point<surface_dim> unit_point = StaticMappingQ1<surface_dim>::mapping
                .transform_real_to_unit_cell(spatial_cell, point);

              spatial_cell->get_dof_values(topography_spatial, topo_dof_values);
              double topo = 0;
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                topo += topo_spatial_fe.shape_value(i, unit_point) * topo_dof_values[i];

              topography_material[z_cell->vertex_dof_index(v, 0)] = topo;
            }
        }
      }


      template <int dim>
      void
      ConvectionDiffusion<dim>::interpolate_coordinates_onto_interface()
      {
        if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) != 0)
          return;

        std::vector<bool> vertex_visited(interface_mesh.n_vertices(), false);
        for (const auto &interface_cell : interface_mesh.active_cell_iterators())
        {
          typename Triangulation<surface_dim>::cell_iterator
          surface_cell = surface_mesh.create_cell_iterator(interface_cell->id());

          for (const unsigned int v : interface_cell->vertex_indices())
            if (!vertex_visited[interface_cell->vertex_index(v)])
            {
              vertex_visited[interface_cell->vertex_index(v)] = true;

              Point<dim> coordinates;
              if (parameters.smooth_surface_mesh)
              {
                // If the surface mesh has been smoothed, then x(y)-coordinates
                // and z-coordinate are stored in vector surface_coordinates and
                // vector topography_material, respectively.
                typename DoFHandler<surface_dim>::cell_iterator
                xy_cell(*surface_cell, &coord_dof_handler);
                while (xy_cell->has_children())
                  xy_cell = xy_cell->child(v);

                typename DoFHandler<surface_dim>::active_cell_iterator
                z_cell(&surface_mesh, xy_cell->level(), xy_cell->index(), 
                       &topo_material_dof_handler);

                for (unsigned int d = 0; d < surface_dim; ++d)
                  coordinates[d] = surface_coordinates[xy_cell->vertex_dof_index(v, d)];

                coordinates[surface_dim] = topography_material[z_cell->vertex_dof_index(v, 0)];

                interface_coordinates[interface_cell->vertex_index(v)] = coordinates;
              }
              else
              {
                // If the surface mesh has not been smoothed, then x(y)-coordinates
                // are identical to the vertex coordinates of the interface mesh, 
                // and z-coordinate is stored in vector topography_spatial.
                typename DoFHandler<surface_dim>::cell_iterator
                z_cell(*surface_cell, &topo_spatial_dof_handler);
                while (z_cell->has_children())
                  z_cell = z_cell->child(v);

                for (unsigned int d = 0; d < surface_dim; ++d)
                  coordinates[d] = interface_cell->vertex(v)(d);

                coordinates[surface_dim] = topography_spatial[z_cell->vertex_dof_index(v, 0)];
              }

              interface_coordinates[interface_cell->vertex_index(v)] = coordinates;
            }
        }
      }


      template <int dim>
      void
      ConvectionDiffusion<dim>::scatter_data_to_surface_owners()
      {
        if (surface_mpi_comm == MPI_COMM_NULL)
          return;

        std::vector<std::vector<char>> send_data_per_process;
        std::vector<void *> send_data_it;
        std::vector<int> send_data_length;

        // The data will be reorganized into a single array.
        std::vector<char>send_data;
        std::vector<int> offsets;

        // For each vertex of the interface mesh, we need to transfer
        // its index and coordinates.
        const unsigned int vertex_size = sizeof(unsigned int) + sizeof(double) * dim;

        // Prepare the data to be sent on root process.
        if (Utilities::MPI::this_mpi_process(surface_mpi_comm) == 0)
        {
          const unsigned int fs_comm_size = Utilities::MPI::n_mpi_processes(surface_mpi_comm);
          send_data_per_process.resize(fs_comm_size);
          send_data_it.resize(fs_comm_size);
          send_data_length.resize(fs_comm_size);

          // For each process in surface_mpi_comm, we need to know whether a vertex
          // has been visited in order to avoid data duplication.
          std::vector<std::vector<bool>> vertex_visited(fs_comm_size);

          for (unsigned int i = 0; i < fs_comm_size; ++i)
          {
            send_data_per_process[i].resize(vertex_size * interface_mesh.n_vertices());
            send_data_it[i] = static_cast<void *>(&send_data_per_process[i].front());
            vertex_visited[i].resize(interface_mesh.n_vertices(), false);
          }

          for (const auto &cell : interface_mesh.active_cell_iterators())
          {
            const unsigned int surface_mpi_rank = interface_to_volume_mapping[cell].surface_mpi_rank;
            for (const unsigned int v : cell->vertex_indices())
            {
              const unsigned int vertex_index = cell->vertex_index(v);
              if (!vertex_visited[surface_mpi_rank][vertex_index])
              {
                vertex_visited[surface_mpi_rank][vertex_index] = true;

                // Write the vertex index.
                unsigned int *index_data = static_cast<unsigned int *>(send_data_it[surface_mpi_rank]);
                *index_data++ = vertex_index;

                // Write the coordinates of this vertex.
                double *coord_data = reinterpret_cast<double *>(index_data);
                for (unsigned int d=0; d<dim; ++d)
                  *coord_data++ = interface_coordinates[vertex_index][d];

                send_data_it[surface_mpi_rank] = static_cast<void *>(coord_data);
              }
            }
          }

          for (unsigned int i = 0; i < fs_comm_size; ++i)
          {
            send_data_length[i] = static_cast<char *>(send_data_it[i]) - &send_data_per_process[i][0];
            send_data.insert(send_data.end(),
                             &send_data_per_process[i][0],
                             static_cast<char *>(send_data_it[i]));
          }

          // Calculate the offset of each segment of data to be sent.
          offsets.push_back(0);
          for (unsigned int i = 1; i < fs_comm_size; ++i)
            offsets.push_back (offsets[i-1] + send_data_length[i-1]);
        }

        // Scatter the data length to each process in surface_mpi_comm.
        int recv_data_length;
        int ierr = MPI_Scatter (&send_data_length[0], 1, MPI_INT, 
                                &recv_data_length, 1, MPI_INT, 
                                0, surface_mpi_comm);
        AssertThrowMPI(ierr);

        Assert(recv_data_length % vertex_size == 0, ExcInternalError());

        // Scatter the coordinate data to each process in surface_mpi_comm.
        std::vector<char> recv_data(recv_data_length);
        ierr = MPI_Scatterv(&send_data[0], &send_data_length[0], &offsets[0], MPI_CHAR,
                            &recv_data[0], recv_data_length, MPI_CHAR,
                            0, surface_mpi_comm);
        AssertThrowMPI(ierr);

        // If this process is the root, it can rest now since it has
        // already had the complete set of elevation values.
        if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
          return;

        // Read the received data.
        const unsigned int n_vertices = recv_data_length / vertex_size;
        const void *recv_data_it = static_cast<const void *>(&recv_data[0]);

        for (unsigned int i = 0; i < n_vertices; ++i)
        {
   // Read the index of the vertex.
          const unsigned int *index_data = static_cast<const unsigned int *>(recv_data_it);
          const unsigned int vertex_index = *index_data++;

          // Read the coordinate values.
          const double *coord_data = reinterpret_cast<const double *>(index_data);
          for (unsigned int d = 0; d < dim; ++d)
            interface_coordinates[vertex_index][d] = *coord_data++;
          recv_data_it = static_cast<const void *>(coord_data);
        }
      }


      template <int dim>
      void ConvectionDiffusion<dim>::initialize()
      {
        // Create the surface mesh.
        const std::set<types::boundary_id> fs_boundary_ids =
          this->get_mesh_deformation_handler().get_free_surface_boundary_indicators();
        internal::extract_coarse_surface_mesh(this->get_triangulation(),
                                              surface_mesh,
                                              fs_boundary_ids);

        unsigned int max_refinement_level = this->get_parameters().initial_global_refinement +
                                            this->get_parameters().initial_adaptive_refinement;
        surface_mesh.refine_global (max_refinement_level);

        global_Omega_diameter = GridTools::diameter(surface_mesh);

        // Get the boundary ids of surface mesh.
        const std::vector<types::boundary_id> 
        surface_mesh_boundary_ids = surface_mesh.get_boundary_ids();

        // Initialize finite elements corresponding to velocity field.
        interface_vel_fe = std::make_unique<FESystem<surface_dim>>(
          FE_Q<surface_dim>(this->get_parameters().displacement_degree), dim);
        surface_vel_fe = std::make_unique<FESystem<surface_dim>>(
          FE_Q<surface_dim>(this->get_parameters().displacement_degree), dim);

        // Initialize the convection-diffusion system.
        {
          // FE space of velocity field.
          surface_vel_dof_handler.distribute_dofs(*surface_vel_fe);
          surface_velocities.reinit(surface_vel_dof_handler.n_dofs());

          // FE space of topography in spatial coordinates.
          topo_spatial_dof_handler.distribute_dofs(topo_spatial_fe);
          topography_spatial.reinit(topo_spatial_dof_handler.n_dofs());
          old_topography_spatial.reinit(topo_spatial_dof_handler.n_dofs());
          old_old_topography_spatial.reinit(topo_spatial_dof_handler.n_dofs());

          // Make boundary constraints for the convection-diffusion system.
          // TODO: Only the zero velocity boundaries are taken into account.
          // However, there exists a special case: some sort of inhomogeneous 
          // Dirichlet boundary conditions are applied to lateral boundaries, and
          // only the velocity in vertical direction is constrained. The 
          // corresponding boundary conditions for the convection-diffusion system
          // should be determined by the velocity boundary model and may change
          // with time. We ignore this case for now because it is rarely met in 
          // practice.
          convection_diffusion_constraints.clear();

          const BoundaryVelocity::Manager<dim> &bv_manager = this->get_boundary_velocity_manager();
          for (const auto boundary_id : bv_manager.get_zero_boundary_velocity_indicators())
            if (std::find(surface_mesh_boundary_ids.begin(), surface_mesh_boundary_ids.end(), boundary_id) != surface_mesh_boundary_ids.end())
            {
              ScalarFunctionFromFunctionObject<surface_dim> scalar_function_object(
                [&] (const Point<surface_dim> &x) -> double
                {
                  return this->get_initial_topography().value(x);
                });

              VectorTools::interpolate_boundary_values(topo_spatial_dof_handler,
                                                       boundary_id,
                                                       scalar_function_object,
                                                       convection_diffusion_constraints);
            }

          convection_diffusion_constraints.close();

          // the linear system of convection-diffusion equation
          DynamicSparsityPattern dsp(topo_spatial_dof_handler.n_dofs(),
                                     topo_spatial_dof_handler.n_dofs());
          DoFTools::make_sparsity_pattern(topo_spatial_dof_handler, dsp,
                                          convection_diffusion_constraints,
                                          false);

          convection_diffusion_sp.copy_from(dsp);
          convection_diffusion_matrix.reinit(convection_diffusion_sp);
          convection_diffusion_rhs.reinit(topo_spatial_dof_handler.n_dofs());
        }

        // Initialize the Laplace-Beltrami system.
        {
          // FE space of topography in material coordinates.
          topo_material_dof_handler.distribute_dofs(topo_material_fe);
          topography_material.reinit(topo_material_dof_handler.n_dofs());

          // FE space of vertex coordinates.
          coord_dof_handler.distribute_dofs(coord_fe);
          surface_coordinates.reinit(coord_dof_handler.n_dofs());

          // Make boundary constraints for the Laplace-Beltrami system.
          // All boundaries are considered to be no-normal-flux boundaries,
          // i.e., mesh vertices are only allowed to move tangentially
          // at boundary.
          internal::make_no_normal_flux_constraints(coord_dof_handler, 
                                                    surface_mesh_boundary_ids,
                                                    laplace_beltrami_constraints);

          // the linear system of Laplace-Beltrami equation
          DynamicSparsityPattern dsp(coord_dof_handler.n_dofs(),
                                     coord_dof_handler.n_dofs());
          DoFTools::make_sparsity_pattern(coord_dof_handler, dsp,
                                          laplace_beltrami_constraints,
                                          false);

          laplace_beltrami_sp.copy_from(dsp);
          laplace_beltrami_matrix.reinit(laplace_beltrami_sp);
          laplace_beltrami_rhs.reinit(coord_dof_handler.n_dofs());
        }

        // Initialize the spatial coordinates (including the vertex coordinates of
        // surface mesh and the topography).
        Quadrature<surface_dim> quadrature(topo_spatial_fe.get_unit_support_points());
        FEValues<surface_dim> fe_values(topo_spatial_fe, quadrature, update_quadrature_points);

        Vector<double> xy_dof_values(coord_fe.dofs_per_cell);
        Vector<double> z_dof_values(topo_spatial_fe.dofs_per_cell);

        typename DoFHandler<surface_dim>::active_cell_iterator
        xy_cell = coord_dof_handler.begin_active(),
        xy_endc = coord_dof_handler.end();
        typename DoFHandler<surface_dim>::active_cell_iterator
        z_cell = topo_spatial_dof_handler.begin_active();
        for (; xy_cell != xy_endc; ++xy_cell, ++z_cell)
        {
          fe_values.reinit(z_cell);
          const std::vector<Point<surface_dim>> &support_points = fe_values.get_quadrature_points();
          for (unsigned int i = 0; i < z_dof_values.size(); ++i)
            z_dof_values[i] = this->get_initial_topography().value(support_points[i]);

          for (const unsigned int v : xy_cell->vertex_indices())
            for (unsigned int d = 0; d < surface_dim; ++d)
            {
              const unsigned int index = coord_fe.component_to_system_index(d, v);
              xy_dof_values[index] = xy_cell->vertex(v)(d);
            }

          xy_cell->set_dof_values(xy_dof_values, surface_coordinates);
          z_cell->set_dof_values(z_dof_values, topography_spatial);
        }

        old_topography_spatial = topography_spatial;
        old_old_topography_spatial = topography_spatial;
      }


      template <int dim>
      void ConvectionDiffusion<dim>::update()
      {
        // The partition may have changed after mesh refinement, so we have to
        // rebuild the sub-communicator composed by surface owners. The complexity
        // of coarsening the interface mesh to the base level and then refining it
        // to the current state is not lower than reconstructing it from beginning
        // (we have to do the MPI communications anyway), so we simply clear the 
        // interface mesh and reconstruct it.
        // TODO: consider the situation that the volume mesh is not refined every
        // time step.
        create_interface_mesh();

        // Gather the velocity data to the root process and distribute them onto 
        // the interface mesh.
        gather_data_to_root_process();

        // The interface mesh is inconsistent with the surface mesh (the former
        // deforms with the volume mesh, while the later stays still), so we 
        // have to execute a projection to transfer velocity data onto the 
        // surface mesh.
        project_velocities_onto_surface();

        // Solve for the topographic evolution, which is considered to be a 
        // convection-diffusion problem.
        compute_topographic_evolution();
        
        // Execute mesh smoothing to the surface mesh if requested by the user.
        if (parameters.smooth_surface_mesh)
          execute_mesh_smoothing();

        // Interpolate the vertex coordinates onto the surface mesh.
        interpolate_coordinates_onto_interface();

        // Scatter the vertex coordinates to the surface owners.
        scatter_data_to_surface_owners();
      }


      template <int dim>
      void
      ConvectionDiffusion<dim>::
      make_boundary_constraints(const TrilinosWrappers::MPI::Vector &/*mesh_displacements*/,
                                const DoFHandler<dim>               &mesh_deformation_dof_handler,
                                AffineConstraints<double>           &mesh_deformation_constraints,
                                const std::set<types::boundary_id>  &fs_boundary_ids) const
      {
        const IndexSet &mesh_locally_owned = mesh_deformation_dof_handler.locally_owned_dofs();
        IndexSet mesh_locally_relevant;
        DoFTools::extract_locally_relevant_dofs(mesh_deformation_dof_handler,
                                                mesh_locally_relevant);

        TrilinosWrappers::MPI::Vector solution, distributed_solution;
        solution.reinit(mesh_locally_owned, mesh_locally_relevant, this->get_mpi_communicator());
        distributed_solution.reinit(mesh_locally_owned, this->get_mpi_communicator());

        const Triangulation<dim> &volume_mesh = mesh_deformation_dof_handler.get_triangulation();
        const FiniteElement<dim> &mesh_deformation_fe = mesh_deformation_dof_handler.get_fe();
        std::vector<types::global_dof_index> cell_dof_indices(mesh_deformation_fe.dofs_per_cell);

        std::vector<bool> vertex_visited(interface_mesh.n_vertices(), false);

        // Loop over the interface mesh and compute the mesh displacements for each vertex.
        for (auto p = volume_to_interface_mapping.begin(); p != volume_to_interface_mapping.end(); ++p)
        {
          const typename Triangulation<dim>::cell_iterator 
          volume_cell = volume_mesh.create_cell_iterator(p->first.adjacent_cell);
          Assert(volume_cell->is_locally_owned(), ExcInternalError());
          
          const unsigned int face_no = p->first.face_no;
          const typename DoFHandler<dim>::cell_iterator dof_cell(*volume_cell,
                                                                 &mesh_deformation_dof_handler);
          const typename Triangulation<surface_dim>::active_cell_iterator &interface_cell = p->second;

          dof_cell->get_dof_indices(cell_dof_indices);
          for (const unsigned int v : interface_cell->vertex_indices())
          {
            const unsigned int vertex_index = interface_cell->vertex_index(v);
            if (!vertex_visited[vertex_index])
            {
              vertex_visited[vertex_index] = true;

              // Get the initial coordinates of the present vertex using function
              // Mapping<dim>::transform_unit_to_real_cell().
              const Point<dim> unit_cell_vertex =
                GeometryInfo<dim>::unit_cell_vertex(GeometryInfo<dim>::face_to_cell_vertices(
                  face_no, v,
                  volume_cell->face_orientation(face_no),
                  volume_cell->face_flip(face_no),
                  volume_cell->face_rotation(face_no)));

              const Point<dim> initial_coordinates = 
                StaticMappingQ1<dim>::mapping.transform_unit_to_real_cell(volume_cell, 
                                                                          unit_cell_vertex);

              for (unsigned int d = 0; d < dim; ++d)
              {
                // Hard coding: the d-th DoF of the v-th vertex is the 
                // (v * dim + d)-th DoF of the face.
                // TODO: the newest version of deal.II has replaced the three boolean arguments
                // of function face_to_cell_index() that determine the orientation of the face
                // by a char type argument that combines the booleans, but I haven't figured out
                // how to get the combined orientation. For now, I simply pass the default value
                // to the function, but it may be wrong for some geometry models other than box.
                const types::global_dof_index dof_index =
                  cell_dof_indices[mesh_deformation_fe.face_to_cell_index(
                    v * dim + d, face_no,
#if DEAL_II_VERSION_GTE(9,6,0)
                    volume_cell->combined_face_orientation(face_no)
#else
                    volume_cell->face_orientation(face_no),
                    volume_cell->face_flip(face_no),
                    volume_cell->face_rotation(face_no)
#endif
                  )];

                if (d < surface_dim)
                {
                  distributed_solution[dof_index] = 
                    interface_coordinates[vertex_index](d) - initial_coordinates(d);
                }
                else
                {
                  Point<surface_dim> position;
                  for (unsigned int j = 0; j < surface_dim; ++j)
                    position(j) = interface_coordinates[vertex_index](j);

                  distributed_solution[dof_index] =
                    interface_coordinates[vertex_index](d) - this->get_initial_topography().value(position);
                }
              }
            }
          }
        }

        solution = distributed_solution;

        // Finally, make constraints for the free surface.
        const IndexSet constrained_dofs =
        DoFTools::extract_boundary_dofs(mesh_deformation_dof_handler,
                                        ComponentMask(dim, true),
                                        fs_boundary_ids);

        for (unsigned int i = 0; i < constrained_dofs.n_elements(); ++i)
        {
          types::global_dof_index index = constrained_dofs.nth_index_in_set(i);
          if (mesh_deformation_constraints.can_store_line(index))
            if (mesh_deformation_constraints.is_constrained(index) == false)
            {
              mesh_deformation_constraints.add_line(index);
              mesh_deformation_constraints.set_inhomogeneity(index, solution[index]);
            }
        }

        mesh_deformation_constraints.close();
      }


      template <int dim>
      void
      ConvectionDiffusion<dim>::declare_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Mesh deformation");
        {
          prm.enter_subsection("Free surface");
          {
            prm.enter_subsection("Convection diffusion");
            {
              prm.declare_entry ("CFL number", "1.0",
                                 Patterns::Double(0),
                                 "");

              prm.declare_entry ("Minimum substep number", "2",
                                 Patterns::Integer(1),
                                 "");

              prm.declare_entry ("Diffusion constant", "0.0",
                                 Patterns::Double(0),
                                 "");

              prm.declare_entry ("cR", "0.33",
                                 Patterns::Double(0),
                                 "");

              prm.declare_entry ("alpha", "2",
                                 Patterns::Integer(1,2),
                                 "");

              prm.declare_entry ("beta", "0.052",
                                 Patterns::Double(0),
                                 "");

              prm.declare_entry ("Smooth surface mesh", "true",
                                 Patterns::Bool(),
                                 "");

              prm.declare_entry ("Gradient descent iterations", "100",
                                 Patterns::Integer(1),
                                 "");

              prm.declare_entry ("Gradient descent step size", "0.01",
                                 Patterns::Double(0),
                                 "");

              prm.declare_entry ("Upper limit of surface slope", "10",
                                 Patterns::Double(0),
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
      ConvectionDiffusion<dim>::parse_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Mesh deformation");
        {
          prm.enter_subsection("Free surface");
          {
            prm.enter_subsection("Convection diffusion");
            {
              parameters.CFL_number                  = prm.get_double("CFL number");
              parameters.diffusion_constant          = prm.get_double("Diffusion constant");
              parameters.minimum_substep_number      = prm.get_integer("Minimum substep number");

              parameters.stabilization_c_R           = prm.get_double("cR");
              parameters.stabilization_alpha         = prm.get_double("alpha");
              parameters.stabilization_beta          = prm.get_double("beta");

              parameters.smooth_surface_mesh         = prm.get_bool("Smooth surface mesh");
              parameters.gradient_descent_iterations = prm.get_integer("Gradient descent iterations");
              parameters.gradient_descent_step_size  = prm.get_double("Gradient descent step size");
              parameters.slope_limit                 = prm.get_double("Upper limit of surface slope");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();

        if (this->convert_output_to_years())
          parameters.diffusion_constant /= year_in_seconds;
      }
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace MeshDeformation
  {
    namespace FreeSurface
    {
      ELASPECT_REGISTER_FREE_SURFACE_MODEL(ConvectionDiffusion,
                                       "convection diffusion",
                                       "")
    }
  }
}
