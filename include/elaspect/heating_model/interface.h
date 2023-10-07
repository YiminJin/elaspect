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


#ifndef _elaspect_heating_model_interface_h
#define _elaspect_heating_model_interface_h

#include <elaspect/plugins.h>
#include <elaspect/simulator_access.h>
#include <elaspect/material_model/io_interface.h>

namespace elaspect
{
  namespace HeatingModel
  {
    using namespace dealii;

    struct HeatingModelOutputs
    {
      HeatingModelOutputs(const unsigned int n_points);

      void reset();

      std::vector<double> heating_source_terms;
      
      std::vector<double> lhs_latent_heat_terms;
    };

    template <int dim>
    class Interface
    {
      public:
        virtual ~Interface() = default;

        virtual void initialize();

        virtual void update();

        virtual
        void
        evaluate(const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                 const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                 HeatingModel::HeatingModelOutputs &heating_model_outputs) const = 0;

        static
        void
        declare_parameters(ParameterHandler &prm);

        virtual
        void
        parse_parameters(ParameterHandler &prm);
    };


    template <int dim>
    class Manager : public SimulatorAccess<dim>
    {
      public:
        void update();

        void
        evaluate(const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                 const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                 HeatingModel::HeatingModelOutputs &heating_model_outputs) const;

        static
        void
        declare_parameters(ParameterHandler &prm);

        void
        parse_parameters(ParameterHandler &prm);

        static
        void
        register_heating_model(const std::string &name,
                               const std::string &description,
                               void (*declare_parameters_function) (ParameterHandler &),
                               Interface<dim> *(*factory_function) ());

        const std::vector<std::string> &
        get_active_heating_model_names() const;

        const std::list<std::unique_ptr<Interface<dim>>> &
        get_active_heating_models() const;

      private:
        std::list<std::unique_ptr<Interface<dim>>> heating_model_objects;

        std::vector<std::string> model_names;
    };


    template <int dim>
    std::string
    get_valid_model_names_pattern ();


#define ELASPECT_REGISTER_HEATING_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ELASPECT_REGISTER_HEATING_MODEL_ ## classname \
  { \
    elaspect::internal::Plugins::RegisterHelper<elaspect::HeatingModel::Interface<2>,classname<2>> \
    dummy_ ## classname ## _2d (&elaspect::HeatingModel::Manager<2>::register_heating_model, \
                                name, description); \
    elaspect::internal::Plugins::RegisterHelper<elaspect::HeatingModel::Interface<3>,classname<3>> \
    dummy_ ## classname ## _3d (&elaspect::HeatingModel::Manager<3>::register_heating_model, \
                                name, description); \
  }
  }
}

#endif
