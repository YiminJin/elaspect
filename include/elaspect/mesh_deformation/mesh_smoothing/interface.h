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


#ifndef _elaspect_mesh_deformation_mesh_smoothing_interface_h
#define _elaspect_mesh_deformation_mesh_smoothing_interface_h

#include <elaspect/global.h>
#include <elaspect/plugins.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/trilinos_vector.h>

namespace elaspect
{
  namespace MeshDeformation
  {
    namespace MeshSmoothing
    {
      using namespace dealii;

      template <int dim>
      class Interface
      {
        public:
          virtual ~Interface() = default;

          virtual void initialize();

          virtual void update();

          virtual
          void
          compute_mesh_displacements(const DoFHandler<dim>           &dof_handler,
                                     const AffineConstraints<double> &constraints,
                                     TrilinosWrappers::MPI::Vector   &displacements) const = 0;

          static
          void
          declare_parameters(ParameterHandler &prm);

          virtual
          void 
          parse_parameters(ParameterHandler &prm);
      };

      template <int dim>
      void
      register_mesh_smoothing_model(const std::string &name,
                                    const std::string &description,
                                    void (*declare_parameters_function) (ParameterHandler &),
                                    Interface<dim> *(*factory_function) ());

      template <int dim>
      Interface<dim> *
      create_mesh_smoothing_model(ParameterHandler &prm);

      template <int dim>
      void
      declare_parameters(ParameterHandler &prm);


#define ELASPECT_REGISTER_MESH_SMOOTHING_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ELASPECT_REGISTER_MESH_SMOOTHING_MODEL_ ## classname \
  { \
    elaspect::internal::Plugins::RegisterHelper<elaspect::MeshDeformation::MeshSmoothing::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&elaspect::MeshDeformation::MeshSmoothing::register_mesh_smoothing_model<2>, \
                                name, description); \
    elaspect::internal::Plugins::RegisterHelper<elaspect::MeshDeformation::MeshSmoothing::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&elaspect::MeshDeformation::MeshSmoothing::register_mesh_smoothing_model<3>, \
                                name, description); \
  }

    }
  }
}

#endif
