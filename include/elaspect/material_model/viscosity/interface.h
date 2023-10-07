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


#ifndef _elaspect_material_model_viscosity_interface_h
#define _elaspect_material_model_viscosity_interface_h

#include <elaspect/global.h>
#include <elaspect/plugins.h>
#include <elaspect/material_model/io_interface.h>

namespace elaspect
{
  namespace MaterialModel
  {
    template <int dim> struct MaterialModelInputs;

    namespace Viscosity
    {
      using namespace dealii;

      template <int dim>
      class Interface
      {
        public:
          virtual ~Interface ();

          virtual void initialize ();

          virtual void update ();

          virtual
          double
          compute_viscosity (const MaterialModel::MaterialModelInputs<dim> &in,
                             const unsigned int i,
                             const unsigned int field) const = 0;

          virtual
          FieldDependences::Dependence
          get_field_dependences () const;

          static
          void
          declare_parameters (ParameterHandler &prm);

          virtual
          void
          parse_parameters (ParameterHandler &prm);
      };

      template <int dim>
      void
      register_viscosity_model (const std::string &name,
                                const std::string &description,
                                void (*declare_parameters_function) (ParameterHandler &),
                                Interface<dim> *(*factory_function) ());

      template <int dim>
      Interface<dim> *
      create_viscosity_model (ParameterHandler &prm);

      template <int dim>
      std::string
      get_valid_model_names_pattern ();

      template <int dim>
      void
      declare_parameters (ParameterHandler &prm);

#define ELASPECT_REGISTER_VISCOSITY_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ELASPECT_REGISTER_VISCOSITY_MODEL_ ## classname \
  { \
    elaspect::internal::Plugins::RegisterHelper<elaspect::MaterialModel::Viscosity::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&elaspect::MaterialModel::Viscosity::register_viscosity_model<2>, \
                                name, description); \
    elaspect::internal::Plugins::RegisterHelper<elaspect::MaterialModel::Viscosity::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&elaspect::MaterialModel::Viscosity::register_viscosity_model<3>, \
                                name, description); \
  }
    }
  }
}


#endif
