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


#ifndef _elaspect_gravity_model_h
#define _elaspect_gravity_model_h

#include <elaspect/plugins.h>

#include <deal.II/base/point.h>

namespace elaspect
{
  namespace GravityModel
  {
    using namespace dealii;

    template <int dim>
    class Interface
    {
      public:
        virtual ~Interface ();

        virtual void initialize ();

        virtual void update ();

        virtual Tensor<1,dim> gravity_vector (const Point<dim> &position) const = 0;

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);
    };

    template <int dim>
    void
    register_gravity_model (const std::string &name,
                            const std::string &description,
                            void (*declare_parameters_function) (ParameterHandler &),
                            Interface<dim> *(*factory_function) ());

    template <int dim>
    Interface<dim> *
    create_gravity_model (ParameterHandler &prm);

    template <int dim>
    void
    declare_parameters (ParameterHandler &prm);


#define ELASPECT_REGISTER_GRAVITY_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ELASPECT_REGISTER_GRAVITY_MODEL_ ## classname \
  { \
    elaspect::internal::Plugins::RegisterHelper<elaspect::GravityModel::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&elaspect::GravityModel::register_gravity_model<2>, \
                                name, description); \
    elaspect::internal::Plugins::RegisterHelper<elaspect::GravityModel::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&elaspect::GravityModel::register_gravity_model<3>, \
                                name, description); \
  }
  }
}

#endif
