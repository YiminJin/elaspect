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


#ifndef _elaspect_boundary_heat_flux_interface_h
#define _elaspect_boundary_heat_flux_interface_h

#include <elaspect/global.h>
#include <elaspect/plugins.h>

namespace elaspect
{
  namespace BoundaryHeatFlux
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
        Tensor<1,dim>
        heat_flux (const types::boundary_id boundary_id,
                   const Point<dim> &position,
                   const Tensor<1,dim> &normal_vector) const = 0;

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);
    };

    template <int dim>
    void
    register_boundary_heat_flux (const std::string &name,
                                 const std::string &description,
                                 void (*declare_parameters_function) (ParameterHandler &),
                                 Interface<dim> *(*factory_function) ());

    template <int dim>
    Interface<dim> *
    create_boundary_heat_flux (ParameterHandler &prm);

    template <int dim>
    void
    declare_parameters (ParameterHandler &prm);

#define ELASPECT_REGISTER_BOUNDARY_HEAT_FLUX_MODEL(classname, name, description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ELASPECT_REGISTER_BOUNDARY_HEAT_FLUX_MODEL_ ## classname \
  { \
    elaspect::internal::Plugins::RegisterHelper<elaspect::BoundaryHeatFlux::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&elaspect::BoundaryHeatFlux::register_boundary_heat_flux<2>, \
                                name, description); \
    elaspect::internal::Plugins::RegisterHelper<elaspect::BoundaryHeatFlux::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&elaspect::BoundaryHeatFlux::register_boundary_heat_flux<3>, \
                                name, description); \
  }
  }
}

#endif
