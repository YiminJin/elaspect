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


#ifndef _elaspect_postprocess_visualization_h
#define _elaspect_postprocess_visualization_h

#include <elaspect/postprocess/interface.h>
#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      class Interface
      {
        public:
          virtual
          ~Interface ();
          
          virtual void initialize ();

          virtual void update ();

          static
          void
          declare_parameters (ParameterHandler &prm);

          virtual
          void
          parse_parameters (ParameterHandler &prm);

          virtual
          std::list<std::string>
          required_other_postprocessors () const;
      };


      template <int dim>
      class CellDataVectorCreator : public Interface<dim>
      {
        public:
          ~CellDataVectorCreator () override = default;

          virtual
          std::vector<std::pair<std::string, Vector<float> *> >
          execute () const = 0;
      };


      template <int dim>
      class SurfaceOnlyVisualization
      {
        public:
          /**
           * Destructor. Made `virtual` to ensure that it is possible to
           * test whether a derived class is derived from this class via
           * a `dynamic_cast`.
           */
          virtual
          ~SurfaceOnlyVisualization () = default;
      };
    }


    template <int dim>
    class Visualization : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        Visualization ();

        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

        void 
        update () override;

        std::list<std::string>
        required_other_postprocessors () const override;

        static
        void
        register_visualization_postprocessor (const std::string &name,
                                              const std::string &description,
                                              void (*declare_parameters_function) (ParameterHandler &),
                                              VisualizationPostprocessors::Interface<dim> *(*factory_function) ());

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        double output_interval;

        double last_output_time;

        unsigned int maximum_timesteps_between_outputs;

        unsigned int last_output_timestep;

        unsigned int output_file_number;

        unsigned int group_files;

        std::string temporary_output_location;

        bool interpolate_output;

        bool write_higher_order_output;

        bool write_in_background_thread;

        void set_last_output_time (const double current_time);

        struct OutputHistory
        {
          OutputHistory ();

          ~OutputHistory ();

          std::string last_mesh_file_name;

          std::vector<std::pair<double,std::string> > times_and_pvtu_names;

          std::vector<std::vector<std::string> > output_file_names_by_timestep;

          std::thread background_thread;
        };

        OutputHistory cell_output_history;

        OutputHistory face_output_history;

        static
        void writer (const std::string filename,
                     const std::string temporary_filename,
                     const std::string *file_contents);

        std::list<std::unique_ptr<VisualizationPostprocessors::Interface<dim> > > postprocessors;

        template <typename DataOutType>
        void write_master_files (const DataOutType &data_out,
                                 const std::string &solution_file_prefix,
                                 const std::vector<std::string> &filenames,
                                 OutputHistory                  &output_history) const;

        template <typename DataOutType>
        std::string write_data_out_data(DataOutType   &data_out,
                                        OutputHistory &output_history) const;
    };


  /**
   * Given a class name, a name, and a description for the parameter file for
   * a postprocessor, register it with the elaspect::Postprocess::Manager class.
   *
   * @ingroup Postprocessing
   */
#define ELASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ELASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR_ ## classname \
  { \
    elaspect::internal::Plugins::RegisterHelper<elaspect::Postprocess::VisualizationPostprocessors::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&elaspect::Postprocess::Visualization<2>::register_visualization_postprocessor, \
                                name, description); \
    elaspect::internal::Plugins::RegisterHelper<elaspect::Postprocess::VisualizationPostprocessors::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&elaspect::Postprocess::Visualization<3>::register_visualization_postprocessor, \
                                name, description); \
  }
  }
}

#endif
