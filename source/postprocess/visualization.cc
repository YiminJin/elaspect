#include <elaspect/postprocess/visualization.h>
#include <elaspect/geometry_model/interface.h>
#include <elaspect/utilities.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_faces.h>

#include <boost/lexical_cast.hpp>

namespace elaspect
{
  namespace Postprocess
  {
    namespace internal
    {
      template <int dim>
      class BaseVariablePostprocessor : public DataPostprocessor<dim>, public SimulatorAccess<dim>
      {
        public:
          void
          evaluate_vector_field (const DataPostprocessorInputs::Vector<dim> &input_data,
                                 std::vector<Vector<double>> &computed_quantities) const override
          {
            const unsigned int n_q_points = input_data.solution_values.size();
            const double one_over_dt = 1. / ( this->get_timestep() * 
                                              ( this->convert_output_to_years() ? 
                                                year_in_seconds : 1. ) );

            for (unsigned int q = 0; q < n_q_points; ++q)
              for (unsigned int i = 0; i < computed_quantities[q].size(); ++i)
              {
                if (this->introspection().component_masks.displacement[i])
                  computed_quantities[q][i] = input_data.solution_values[q][i] * one_over_dt;
                else
                  computed_quantities[q][i] = input_data.solution_values[q][i];
              }
          }


          std::vector<std::string> get_names () const override
          {
            std::vector<std::string> solution_names (dim, "velocity");
            solution_names.emplace_back("temperature");

            for (unsigned int c = 0; c < this->n_compositional_fields(); ++c)
              solution_names.emplace_back(this->introspection().name_for_compositional_index(c));

            solution_names.emplace_back("stress_xx");
            solution_names.emplace_back("stress_yy");
            if (dim == 2)
              solution_names.emplace_back("stress_xy");
            else
            {
              solution_names.emplace_back("stress_zz");
              solution_names.emplace_back("stress_xy");
              solution_names.emplace_back("stress_xz");
              solution_names.emplace_back("stress_yz");
            }

            if (this->constitutive_relation() & ConstitutiveRelation::plasticity)
              solution_names.emplace_back("plastic_strain");

            if (this->constitutive_relation() & ConstitutiveRelation::pore_fluid)
              solution_names.emplace_back("fluid_pressure");

            return solution_names;
          }


          std::vector<DataComponentInterpretation::DataComponentInterpretation>
          get_data_component_interpretation () const override
          {
            // displacement
            std::vector<DataComponentInterpretation::DataComponentInterpretation>
            interpretation(dim, DataComponentInterpretation::component_is_part_of_vector);

            // temperature
            interpretation.push_back(DataComponentInterpretation::component_is_scalar);

            // compositional fields
            for (unsigned int c = 0; c < this->n_compositional_fields(); ++c)
              interpretation.push_back(DataComponentInterpretation::component_is_scalar);

            // stress
            for (unsigned int i = 0; i < SymmetricTensor<2,dim>::n_independent_components; ++i)
              interpretation.push_back(DataComponentInterpretation::component_is_scalar);

            if (this->constitutive_relation() & ConstitutiveRelation::plasticity)
              // plastic strain
              interpretation.push_back(DataComponentInterpretation::component_is_scalar);

            if (this->constitutive_relation() & ConstitutiveRelation::pore_fluid)
              // fluid pressure
              interpretation.push_back(DataComponentInterpretation::component_is_scalar);

            return interpretation;
          }

          UpdateFlags get_needed_update_flags () const override
          {
            return update_values;
          }
      };
    }


    namespace VisualizationPostprocessors
    {
      template <int dim>
      Interface<dim>::~Interface ()
      {}


      template <int dim>
      void
      Interface<dim>::initialize ()
      {}


      template <int dim>
      void
      Interface<dim>::update ()
      {}


      template <int dim>
      void
      Interface<dim>::declare_parameters (ParameterHandler &)
      {}



      template <int dim>
      void
      Interface<dim>::parse_parameters (ParameterHandler &)
      {}



      template <int dim>
      std::list<std::string>
      Interface<dim>::required_other_postprocessors () const
      {
        return std::list<std::string>();
      }
    }


    template <int dim>
    Visualization<dim>::OutputHistory::OutputHistory ()
    {}

    
    template <int dim>
    Visualization<dim>::OutputHistory::~OutputHistory ()
    {
      // Make sure that any thread that may still be running in the background,
      // writing data, finishes
      if (background_thread.joinable())
        background_thread.join ();
    }


    template <int dim>
    Visualization<dim>::Visualization ()
      : 
      // the following value is later read from the input file
      output_interval (0),
      // initialize this to a nonsensical value; set it to the actual time
      // the first time around we get to check it
      last_output_time (std::numeric_limits<double>::quiet_NaN()),
      maximum_timesteps_between_outputs (std::numeric_limits<int>::max()),
      last_output_timestep (numbers::invalid_unsigned_int),
      output_file_number (numbers::invalid_unsigned_int)
    {}


    template <int dim>
    void
    Visualization<dim>::update ()
    {
      //Call the .update() method for each visualization postprocessor.
      for (const auto &p : postprocessors)
        p->update();
    }


    template <int dim>
    template <typename DataOutType>
    void
    Visualization<dim>::
    write_master_files (const DataOutType              &data_out,
                        const std::string              &solution_file_prefix,
                        const std::vector<std::string> &filenames,
                        OutputHistory                  &output_history) const
    {
      static_assert (std::is_same<DataOutType,DataOut<dim>>::value ||
                     std::is_same<DataOutType,DataOutFaces<dim>>::value,
                     "The only allowed template types of this function are "
                     "DataOut and DataOutFaces.");
      const bool is_cell_data_output = std::is_same<DataOutType,DataOut<dim> >::value;

      const double time_in_years_or_seconds = (this->convert_output_to_years() ?
                                               this->get_time() / year_in_seconds :
                                               this->get_time());
      const std::string pvtu_master_filename = (solution_file_prefix + ".pvtu");

      std::ofstream pvtu_master ((this->get_output_directory() + "solution/" +
                                  pvtu_master_filename).c_str());
      data_out.write_pvtu_record (pvtu_master, filenames);

      // now also generate a .pvd file that matches simulation
      // time and corresponding .pvtu record
      if (this->get_parameters().run_postprocessors_on_nonlinear_iterations)
      {
        // in case we output all nonlinear iterations, we only want one
        // entry per time step, so replace the last line with the current iteration
        if (this->get_nonlinear_iteration() == 0)
          output_history.times_and_pvtu_names.emplace_back(time_in_years_or_seconds, "solution/"+pvtu_master_filename);
        else
          output_history.times_and_pvtu_names.back() = (std::make_pair
                                                        (time_in_years_or_seconds, "solution/"+pvtu_master_filename));
      } 
      else
        output_history.times_and_pvtu_names.emplace_back(time_in_years_or_seconds, 
                                                         "solution/"+pvtu_master_filename);

      const std::string pvd_master_filename = (this->get_output_directory() + 
                                               (is_cell_data_output ? "solution.pvd" : "solution_surface.pvd"));
      std::ofstream pvd_master (pvd_master_filename.c_str());

      DataOutBase::write_pvd_record (pvd_master, output_history.times_and_pvtu_names);

      // finally, do the same for Visit via the .visit file for this
      // time step, as well as for all time steps together
      const std::string visit_master_filename = (this->get_output_directory()
                                                 + "solution/"
                                                 + solution_file_prefix
                                                 + ".visit");
      std::ofstream visit_master (visit_master_filename.c_str());

      DataOutBase::write_visit_record (visit_master, filenames);

      {
        // the global .visit file needs the relative path because it sits a
        // directory above
        std::vector<std::string> filenames_with_path;
        filenames_with_path.reserve(filenames.size());
        for (const auto &filename : filenames)
          filenames_with_path.emplace_back("solution/" + filename);

        if (this->get_parameters().run_postprocessors_on_nonlinear_iterations)
        {
          // in case we output all nonlinear iterations, we only want one
          // entry per time step, so replace the last line with the current iteration
          if (this->get_nonlinear_iteration() == 0)
            output_history.output_file_names_by_timestep.push_back (filenames_with_path);
          else
            output_history.output_file_names_by_timestep.back() = filenames_with_path;
        }
        else
          output_history.output_file_names_by_timestep.push_back (filenames_with_path);
      }

      std::ofstream global_visit_master ((this->get_output_directory() + 
                                          (is_cell_data_output ? "solution.visit" : "solution_surface.visit")).c_str());

      std::vector<std::pair<double, std::vector<std::string> > > times_and_output_file_names;
      for (unsigned int timestep=0; timestep<output_history.times_and_pvtu_names.size(); ++timestep)
        times_and_output_file_names.push_back(std::make_pair(output_history.times_and_pvtu_names[timestep].first,
                                                             output_history.output_file_names_by_timestep[timestep]));
      DataOutBase::write_visit_record (global_visit_master, times_and_output_file_names);
    }


    template <int dim>
    template <typename DataOutType>
    std::string
    Visualization<dim>::write_data_out_data (DataOutType   &data_out,
                                             OutputHistory &output_history) const
    {
      static_assert (std::is_same<DataOutType,DataOut<dim>>::value ||
                     std::is_same<DataOutType,DataOutFaces<dim>>::value,
                     "The only allowed template types of this function are "
                     "DataOut and DataOutFaces.");
      const bool is_cell_data_output = std::is_same<DataOutType,DataOut<dim>>::value;

      const double time_in_years_or_seconds = (this->convert_output_to_years() ?
                                               this->get_time() / year_in_seconds :
                                               this->get_time());

      std::string solution_file_prefix = 
        (is_cell_data_output ? "solution-" : "solution_surface-") 
        + Utilities::int_to_string(output_file_number, 5);
      if (this->get_parameters().run_postprocessors_on_nonlinear_iterations)
        solution_file_prefix.append("." + Utilities::int_to_string(this->get_nonlinear_iteration(), 4));

      // Write master files (.pvtu,.pvd,.visit) on the master process
      const int my_id = Utilities::MPI::this_mpi_process(this->get_mpi_communicator());
      if (my_id == 0)
      {
        std::vector<std::string> filenames;
        const unsigned int n_processes = Utilities::MPI::n_mpi_processes(
                                           this->get_mpi_communicator());
        const unsigned int n_files =
          (group_files == 0) ?
          n_processes : std::min(group_files, n_processes);
        for (unsigned int i = 0; i < n_files; ++i)
          filenames.push_back(
            solution_file_prefix + "."
            + Utilities::int_to_string(i, 4) + ".vtu");
        write_master_files(data_out, solution_file_prefix, filenames, output_history);
      }
      
      const unsigned int n_processes = Utilities::MPI::n_mpi_processes(
                                         this->get_mpi_communicator());
      const unsigned int my_file_id = (
                                        group_files == 0 ? my_id : my_id % group_files);
      const std::string filename = this->get_output_directory() + "solution/"
                                   + solution_file_prefix + "."
                                   + Utilities::int_to_string(my_file_id, 4) + ".vtu";

      // pass time step number and time as metadata into the output file
      DataOutBase::VtkFlags vtk_flags;
      vtk_flags.cycle = this->get_timestep_number();
      vtk_flags.time = time_in_years_or_seconds;

      data_out.set_flags(vtk_flags);

      // Write as many files as processes. For this case we support writing in a
      // background thread and to a temporary location, so we first write everything
      // into a string that is written to disk in a writer function
      if ((group_files == 0) || (group_files >= n_processes))
      {
        // Put the content we want to write into a string object that
        // we can then write in the background
        const std::string *file_contents;
        {
          std::ostringstream tmp;
          data_out.write(tmp,
                         DataOutBase::parse_output_format("vtu"));
          file_contents = new std::string(tmp.str());
        }
        if (write_in_background_thread)
        {
          // Wait for all previous write operations to finish, should
          // any be still active, ...
          if (output_history.background_thread.joinable())
            output_history.background_thread.join();
          // ...then continue with writing our own data.
          output_history.background_thread
            = std::thread([&]()
          {
            writer (filename, temporary_output_location, file_contents);
          });
        }
        else
          writer(filename, temporary_output_location, file_contents);
      }
      else
      {
        // Just write one data file in parallel
        if (group_files == 1)
        {
          data_out.write_vtu_in_parallel(filename.c_str(),
                                         this->get_mpi_communicator());
        }
        else               // Write as many output files as 'group_files' groups
        {
          int color = my_id % group_files;
          MPI_Comm comm;
          int ierr = MPI_Comm_split(this->get_mpi_communicator(), color, my_id, &comm);
          AssertThrowMPI(ierr);
          data_out.write_vtu_in_parallel(filename.c_str(), comm);
          ierr = MPI_Comm_free(&comm);
          AssertThrowMPI(ierr);
        }
      }

      return solution_file_prefix;
    }


    namespace
    {
      // Add new output_data_names to a set output_data_names_set, and check for uniqueness of names.
      // Multiple copies in output_data_names are collapsed into one copy before checking
      // for uniqueness with output_data_names_set, because vector data fields are represented as multiple
      // copies of the same name.
      void add_data_names_to_set(const std::vector<std::string> &output_data_names,
                                 std::set<std::string> &output_data_names_set)
      {
        const std::set<std::string> set_of_names(output_data_names.begin(),output_data_names.end());

        for (const auto &name: set_of_names)
          {
            const auto iterator_and_success = output_data_names_set.insert(name);
            AssertThrow(iterator_and_success.second == true,
                        ExcMessage("The output variable <" + name + "> already exists in the list of output "
                                   "variables. Make sure there is no duplication in the names of visualization output "
                                   "variables, otherwise output files may be corrupted."));
          }
      }
    }


    template <int dim>
    std::pair<std::string,std::string>
    Visualization<dim>::execute (TableHandler &statistics)
    {
      // if this is the first time we get here, set the last output time
      // to the current time - output_interval. this makes sure we
      // always produce data during the first time step
      if (std::isnan(last_output_time))
      {
        last_output_time = this->get_time() - output_interval;
        last_output_timestep = this->get_timestep_number();
      }

      // Return if graphical output is not requested at this time. Do not
      // return in the first timestep, or if the last output was more than
      // output_interval in time ago, or maximum_timesteps_between_outputs in
      // number of timesteps ago.
      // The comparison in number of timesteps is safe from integer overflow for
      // at most 2 billion timesteps , which is not likely to
      // be ever reached (both values are unsigned int,
      // and the default value of maximum_timesteps_between_outputs is
      // set to numeric_limits<int>::max())
      if ((this->get_time() < last_output_time + output_interval)
          && (this->get_timestep_number() < last_output_timestep + maximum_timesteps_between_outputs)
          && (this->get_timestep_number() != 0))
        return std::pair<std::string,std::string>();

      // up the counter of the number of the file by one, but not in
      // the very first output step. if we run postprocessors on all
      // iterations, only increase file number in the first nonlinear iteration
      const bool increase_file_number = (this->get_nonlinear_iteration() == 0) || (!this->get_parameters().run_postprocessors_on_nonlinear_iterations);
      if (output_file_number == numbers::invalid_unsigned_int)
        output_file_number = 0;
      else if (increase_file_number)
        ++output_file_number;

      internal::BaseVariablePostprocessor<dim> base_variables;
      base_variables.initialize_simulator (this->get_simulator());

      // Keep a list of the names of all output variables, to ensure unique names
      std::set<std::string> visualization_field_names;

      // Insert base variable names into set of all output field names
      add_data_names_to_set(base_variables.get_names(), visualization_field_names);

      DataOut<dim> data_out;
      data_out.attach_dof_handler (this->get_dof_handler());
      data_out.add_data_vector (this->get_solution(),
                                base_variables);

      // Also create an object for outputting information that lives on
      // the faces of the mesh. If there are postprocessors derived from
      // the VisualizationPostprocessors::SurfaceOnlyVisualization class, then
      // we will use this object for viz purposes.
      DataOutFaces<dim> data_out_faces;
      data_out_faces.attach_dof_handler (this->get_dof_handler());
      const bool have_face_viz_postprocessors
        = (std::find_if (postprocessors.begin(),
                         postprocessors.end(),
                         [](const std::unique_ptr<VisualizationPostprocessors::Interface<dim> > &p)
      {
        return (dynamic_cast<const VisualizationPostprocessors::SurfaceOnlyVisualization<dim>*>
                (p.get()) != nullptr);
      })
      != postprocessors.end());

      // then for each additional selected output variable
      // add the computed quantity as well. keep a list of
      // pointers to data vectors created by cell data visualization
      // postprocessors that will later be deleted
      std::list<std::unique_ptr<Vector<float> > > cell_data_vectors;
      for (const auto &p : postprocessors)
      {
        try
        {
          // There are two ways of writing visualization postprocessors:
          // - deriving from DataPostprocessor
          // - deriving from DataVectorCreator
          // treat them in turn. In both cases, the information can
          // be output on all cells, or via the faces on the surface
          // only (if the class in question is derived from
          // SurfaceOnlyVisualization), so we will have to switch between
          // the two ways when we send things to the data_out or
          // data_out_faces objects.
          if (const DataPostprocessor<dim> *viz_postprocessor
              = dynamic_cast<const DataPostprocessor<dim>*>(& *p))
            {
              add_data_names_to_set(viz_postprocessor->get_names(), visualization_field_names);

              if (dynamic_cast<const VisualizationPostprocessors::SurfaceOnlyVisualization<dim>*>
                  (& *p) == nullptr)
                data_out.add_data_vector (this->get_solution(),
                                          *viz_postprocessor);
              else
                data_out_faces.add_data_vector (this->get_solution(),
                                                *viz_postprocessor);
            }
          else if (const VisualizationPostprocessors::CellDataVectorCreator<dim> *
                   cell_data_creator
                   = dynamic_cast<const VisualizationPostprocessors::CellDataVectorCreator<dim>*>
                     (& *p))
          {
            // get the data produced here
            const std::vector<std::pair<std::string, Vector<float> *> >
            cell_data_set = cell_data_creator->execute();
            for (const auto &cell_data : cell_data_set)
            {
              Assert (cell_data.second->size() ==
                      this->get_triangulation().n_active_cells(),
                      ExcMessage ("Cell data visualization postprocessors must generate "
                                  "vectors that have as many entries as there are active cells "
                                  "on the current processor."));

              add_data_names_to_set(std::vector<std::string>(1,cell_data.first), visualization_field_names);

              // store the pointer, then attach the vector to the DataOut object
              cell_data_vectors.push_back (std::unique_ptr<Vector<float> >
                                           (cell_data.second));

              if (dynamic_cast<const VisualizationPostprocessors::SurfaceOnlyVisualization<dim>*>
                  (& *p) == nullptr)
                data_out.add_data_vector (*cell_data.second,
                                          cell_data.first,
                                          DataOut<dim>::type_cell_data);
              else
                data_out_faces.add_data_vector (*cell_data.second,
                                                cell_data.first,
                                                DataOutFaces<dim>::type_cell_data);
            }
          }
          else
            // A viz postprocessor not derived from either DataPostprocessor
            // or CellDataVectorCreator? We don't know what to do with that!
            Assert (false,
                    ExcMessage("The visualization system found a visualization "
                               "postprocessor class that is not either derived from "
                               "DataPostprocessor or "
                               "VisualizationPostprocessors::CellDataVectorCreator. "
                               "elASPECT does not know what to do with these kinds of "
                               "classes."));
        }
        // viz postprocessors that throw exceptions usually do not result in
        // anything good because they result in an unwinding of the stack
        // and, if only one processor triggers an exception, the
        // destruction of objects often causes a deadlock. thus, if
        // an exception is generated, catch it, print an error message,
        // and abort the program
        catch (std::exception &exc)
        {
          std::cerr << std::endl << std::endl
                    << "----------------------------------------------------"
                    << std::endl;
          std::cerr << "An exception happened on MPI process <"
                    << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
                    << "> while running the visualization postprocessor <"
                    << typeid(*p).name()
                    << ">: " << std::endl
                    << exc.what() << std::endl
                    << "Aborting!" << std::endl
                    << "----------------------------------------------------"
                    << std::endl;

          // terminate the program!
          MPI_Abort (MPI_COMM_WORLD, 1);
        }
        catch (...)
        {
          std::cerr << std::endl << std::endl
                    << "----------------------------------------------------"
                    << std::endl;
          std::cerr << "An exception happened on MPI process <"
                    << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
                    << "> while running the visualization postprocessor <"
                    << typeid(*p).name()
                    << ">: " << std::endl;
          std::cerr << "Unknown exception!" << std::endl
                    << "Aborting!" << std::endl
                    << "----------------------------------------------------"
                    << std::endl;

          // terminate the program!
          MPI_Abort (MPI_COMM_WORLD, 1);
        }
      }

      // Now build the patches. If selected, increase the output resolution.
      // Giving the mapping ensures that the case with mesh deformation works correctly.
      const unsigned int subdivisions = interpolate_output
                                        ? this->get_displacement_degree()
                                        : 0;

      // Now get everything written for the DataOut case, and record this
      // in the statistics file
      std::string solution_file_prefix;
      {
        data_out.build_patches (this->get_mapping(),
                                subdivisions,
                                DataOut<dim>::no_curved_cells);

        solution_file_prefix
          = write_data_out_data(data_out, cell_output_history);
        statistics.add_value ("Visualization file name",
                              this->get_output_directory()
                              + "solution/"
                              + solution_file_prefix);
      }

      // Then do the same again for the face data case. We won't print the
      // output file name to screen (too much clutter on the screen already)
      // but still put it into the statistics file
      if (have_face_viz_postprocessors)
      {
        data_out_faces.build_patches (this->get_mapping(),
                                      subdivisions);

        const std::string face_solution_file_prefix
          = write_data_out_data(data_out_faces, face_output_history);
        statistics.add_value ("Surface visualization file name",
                              this->get_output_directory()
                              + "solution_surface/"
                              + face_solution_file_prefix);
      }

      // Increment the next time we need output:
      set_last_output_time (this->get_time());
      last_output_timestep = this->get_timestep_number();

      // Return what should be printed to the screen. This is a bit
      // late (the output has already been written, and this probably took
      // a good long while), but it's still good to provide a status
      // update.
      return std::make_pair (std::string ("Writing graphical output:"),
                             this->get_output_directory()
                             + "solution/"
                             + solution_file_prefix);
    }


    template <int dim>
    // We need to pass the arguments by value, as this function can be called on a separate thread:
    void Visualization<dim>::writer (const std::string filename, //NOLINT(performance-unnecessary-value-param)
                                     const std::string temporary_output_location, //NOLINT(performance-unnecessary-value-param)
                                     const std::string *file_contents)
    {
      std::string tmp_filename = filename;
      if (temporary_output_location != "")
      {
        tmp_filename = temporary_output_location + "/elaspect.tmp.XXXXXX";

        // Create the temporary file; get at the actual filename
        // by using a C-style string that mkstemp will then overwrite
        std::vector<char> tmp_filename_x (tmp_filename.size()+1);
        std::strcpy(tmp_filename_x.data(), tmp_filename.c_str());
        const int tmp_file_desc = mkstemp(tmp_filename_x.data());
        tmp_filename = tmp_filename_x.data();

        // If we failed to create the temp file, just write directly to the target file.
        // We also provide a warning about this fact. There are places where
        // this fails *on every node*, so we will get a lot of warning messages
        // into the output; in these cases, just writing multiple pieces to
        // std::cerr will produce an unreadable mass of text; rather, first
        // assemble the error message completely, and then output it atomically
        if (tmp_file_desc == -1)
        {
          const std::string x = ("***** WARNING: could not create temporary file <"
                                 +
                                 tmp_filename
                                 +
                                 ">, will output directly to final location. This may negatively "
                                 "affect performance. (On processor "
                                 + Utilities::int_to_string(Utilities::MPI::this_mpi_process (MPI_COMM_WORLD))
                                 + ".)\n");

          std::cerr << x << std::flush;

          tmp_filename = filename;
        }
        else
          close(tmp_file_desc);
      }

      std::ofstream out(tmp_filename.c_str());

      AssertThrow (out, ExcMessage(std::string("Trying to write to file <") +
                                   filename +
                                   ">, but the file can't be opened!"))

      // now write and then move the tmp file to its final destination
      // if necessary
      out << *file_contents;
      out.close ();

      if (tmp_filename != filename)
      {
        std::string command = std::string("mv ") + tmp_filename + " " + filename;
        int error = system(command.c_str());

        AssertThrow(error == 0,
                    ExcMessage("Could not move " + tmp_filename + " to "
                               + filename + ". On processor "
                               + Utilities::int_to_string(Utilities::MPI::this_mpi_process (MPI_COMM_WORLD)) + "."));
      }

      // destroy the pointer to the data we needed to write
      delete file_contents;
    }


    namespace
    {
      std::tuple
      <void *,
      void *,
      elaspect::internal::Plugins::PluginList<VisualizationPostprocessors::Interface<2> >,
      elaspect::internal::Plugins::PluginList<VisualizationPostprocessors::Interface<3> > > registered_visualization_plugins;
    }


    template <int dim>
    void
    Visualization<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Visualization");
        {
          prm.declare_entry ("Time between graphical output", "1e8",
                             Patterns::Double (0.),
                             "The time interval between each generation of "
                             "graphical output files. A value of zero indicates "
                             "that output should be generated in each time step. "
                             "Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");

          prm.declare_entry ("Time steps between graphical output", boost::lexical_cast<std::string>(std::numeric_limits<int>::max()),
                             Patterns::Integer(0),
                             "The maximum number of time steps between each generation of "
                             "graphical output files.");

          // now also see about the file format we're supposed to write in
          prm.declare_entry ("Output format", "vtu",
                             Patterns::Selection (DataOutBase::get_output_format_names ()),
                             "The file format to be used for graphical output. The list "
                             "of possible output formats that can be given here is documented "
                             "in the appendix of the manual where the current parameter "
                             "is described.");

          prm.declare_entry ("Number of grouped files", "16",
                             Patterns::Integer(0),
                             "VTU file output supports grouping files from several CPUs "
                             "into a given number of files using MPI I/O when writing on a parallel "
                             "filesystem. Select 0 for no grouping. This will disable "
                             "parallel file output and instead write one file per processor. "
                             "A value of 1 will generate one big file containing the whole "
                             "solution, while a larger value will create that many files "
                             "(at most as many as there are MPI ranks).");

          prm.declare_entry ("Write in background thread", "false",
                             Patterns::Bool(),
                             "File operations can potentially take a long time, blocking the "
                             "progress of the rest of the model run. Setting this variable to "
                             "`true' moves this process into a background thread, while the "
                             "rest of the model continues.");

          prm.declare_entry ("Temporary output location", "",
                             Patterns::Anything(),
                             "On large clusters it can be advantageous to first write the "
                             "output to a temporary file on a local file system and later "
                             "move this file to a network file system. If this variable is "
                             "set to a non-empty string it will be interpreted as a "
                             "temporary storage location.");

          prm.declare_entry ("Interpolate output", "true",
                             Patterns::Bool(),
                             "deal.II offers the possibility to linearly interpolate "
                             "output fields of higher order elements to a finer resolution. "
                             "This somewhat compensates the fact that most visualization "
                             "software only offers linear interpolation between grid points "
                             "and therefore the output file is a very coarse representation "
                             "of the actual solution field. Activating this option increases "
                             "the spatial resolution in each dimension by a factor equal "
                             "to the polynomial degree used for the velocity finite element "
                             "(usually 2). In other words, instead of showing one quadrilateral "
                             "or hexahedron in the visualization per cell on which \\elaspect{} "
                             "computes, it shows multiple (for quadratic elements, it will "
                             "describe each cell of the mesh on which we compute as "
                             "$2\\times 2$ or $2\\times 2\\times 2$ cells in 2d and 3d, "
                             "respectively; correspondingly more subdivisions are used if "
                             "you use cubic, quartic, or even higher order elements for the "
                             "velocity)."
                             "\n\n"
                             "The effect of using this option can be seen in the following "
                             "picture showing a variation of the output produced with the "
                             "input files from Section~\\ref{sec:shell-simple-2d}:"
                             "\n\n"
                             "\\begin{center}"
                             "  \\includegraphics[width=0.5\\textwidth]{viz/parameters/build-patches}"
                             "\\end{center}"
                             "Here, the left picture shows one visualization cell per "
                             "computational cell (i.e., the option is switched off), "
                             "and the right picture shows the same simulation with the "
                             "option switched on (which is the default). The images "
                             "show the same data, demonstrating "
                             "that interpolating the solution onto bilinear shape functions as is "
                             "commonly done in visualizing data loses information."
                             "\n\n"
                             "Of course, activating this option also greatly increases the amount of "
                             "data \\elaspect{} will write to disk: approximately by a factor of 4 in 2d, "
                             "and a factor of 8 in 3d, when using quadratic elements for the velocity, "
                             "and correspondingly more for even higher order elements.");

          prm.declare_entry ("Write higher order output", "false",
                             Patterns::Bool(),
                             "deal.II offers the possibility to write vtu files with higher order "
                             "representations of the output data. This means each cell will correctly "
                             "show the higher order representation of the output data instead of the "
                             "linear interpolation between vertices that ParaView and Visit usually show. "
                             "Note that activating this option is safe and recommended, but requires that "
                             "(i) ``Output format'' is set to ``vtu'', (ii) ``Interpolate output'' is "
                             "set to true, (iii) you use a sufficiently new version of Paraview "
                             "or Visit to read the files (Paraview version 5.5 or newer, and Visit version "
                             "to be determined), and (iv) you use deal.II version 9.1.0 or newer. "
                             "\n"
                             "The effect of using this option can be seen in the following "
                             "picture:"
                             "\n\n"
                             "\\begin{center}"
                             "  \\includegraphics[width=0.5\\textwidth]{viz/parameters/higher-order-output}"
                             "\\end{center}"
                             "The top figure shows the plain output without interpolation or higher "
                             "order output. The middle figure shows output that was interpolated as "
                             "discussed for the ``Interpolate output'' option. The bottom panel "
                             "shows higher order output that achieves better accuracy than the "
                             "interpolated output at a lower memory cost.");

          // Finally also construct a string for Patterns::MultipleSelection that
          // contains the names of all registered visualization postprocessors.
          // Also add a number of removed plugins that are now combined in 'material properties'
          // to keep compatibility with input files. These will be filtered out in parse_parameters().
          const std::string pattern_of_names
            = std::get<dim>(registered_visualization_plugins).get_pattern_of_names ()
              + "|density|specific heat|thermal conductivity|thermal diffusivity|thermal expansivity|viscosity";
          prm.declare_entry("List of output variables",
                            "",
                            Patterns::MultipleSelection(pattern_of_names),
                            "A comma separated list of visualization objects that should be run "
                            "whenever writing graphical output. By default, the graphical "
                            "output files will always contain the primary variables velocity, "
                            "pressure, and temperature. However, one frequently wants to also "
                            "visualize derived quantities, such as the thermodynamic phase "
                            "that corresponds to a given temperature-pressure value, or the "
                            "corresponding seismic wave speeds. The visualization objects do "
                            "exactly this: they compute such derived quantities and place them "
                            "into the output file. The current parameter is the place where "
                            "you decide which of these additional output variables you want "
                            "to have in your output file.\n\n"
                            "The following postprocessors are available:\n\n"
                            +
                            std::get<dim>(registered_visualization_plugins).get_description_string());
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // now declare the parameters of each of the registered
      // visualization postprocessors in turn
      std::get<dim>(registered_visualization_plugins).declare_parameters (prm);
    }


    template <int dim>
    void
    Visualization<dim>::parse_parameters (ParameterHandler &prm)
    {
      Assert (std::get<dim>(registered_visualization_plugins).plugins != nullptr,
              ExcMessage ("No postprocessors registered!?"));
      std::vector<std::string> viz_names;

      std::string visualization_subdirectory = this->get_output_directory() + "solution/";
      Utilities::create_directory (visualization_subdirectory,
                                   this->get_mpi_communicator(),
                                   true);

      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Visualization");
        {
          output_interval = prm.get_double ("Time between graphical output");
          if (this->convert_output_to_years())
            output_interval *= year_in_seconds;

          maximum_timesteps_between_outputs = prm.get_integer("Time steps between graphical output");

          if (output_interval > 0.0)
          {
            // since we increase the time indicating when to write the next graphical output
            // every time we execute the visualization postprocessor, there is no good way to
            // figure out when to write graphical output for the nonlinear iterations if we do
            // not want to output every time step
            AssertThrow(this->get_parameters().run_postprocessors_on_nonlinear_iterations == false,
                        ExcMessage("Postprocessing nonlinear iterations is only supported if every time "
                                   "step is visualized, or in other words, if the 'Time between graphical "
                                   "output' in the Visualization postprocessor is set to zero."));
          }

          group_files     = prm.get_integer("Number of grouped files");
          write_in_background_thread = prm.get_bool("Write in background thread");
          temporary_output_location = prm.get("Temporary output location");

          if (temporary_output_location != "")
          {
            // Check if a command-processor is available by calling system() with a
            // null pointer. System is guaranteed to return non-zero if it finds
            // a terminal and zero if there is none (like on the compute nodes of
            // some cluster architectures, e.g. IBM BlueGene/Q)
            AssertThrow(system((char *)nullptr) != 0,
                        ExcMessage("Usage of a temporary storage location is only supported if "
                                   "there is a terminal available to move the files to their final location "
                                   "after writing. The system() command did not succeed in finding such a terminal."));
          }

          interpolate_output = prm.get_bool("Interpolate output");

          // now also see which derived quantities we are to compute
          viz_names = Utilities::split_string_list(prm.get("List of output variables"));
          AssertThrow(Utilities::has_unique_entries(viz_names),
                      ExcMessage("The list of strings for the parameter "
                                 "'Postprocess/Visualization/List of output variables' contains entries more than once. "
                                 "This is not allowed. Please check your parameter file."));

          // see if 'all' was selected (or is part of the list). if so
          // simply replace the list with one that contains all names
          if (std::find (viz_names.begin(),
                         viz_names.end(),
                         "all") != viz_names.end())
          {
            viz_names.clear();
            for (typename std::list<typename elaspect::internal::Plugins::PluginList<VisualizationPostprocessors::Interface<dim> >::PluginInfo>::const_iterator
                 p = std::get<dim>(registered_visualization_plugins).plugins->begin();
                 p != std::get<dim>(registered_visualization_plugins).plugins->end(); ++p)
              viz_names.push_back (std::get<0>(*p));
          }
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // then go through the list, create objects and let them parse
      // their own parameters
      for (unsigned int name=0; name<viz_names.size(); ++name)
      {
        VisualizationPostprocessors::Interface<dim> *
        viz_postprocessor = std::get<dim>(registered_visualization_plugins)
                            .create_plugin (viz_names[name],
                                            "Visualization plugins");

        // make sure that the postprocessor is indeed of type
        // dealii::DataPostprocessor or of type
        // VisualizationPostprocessors::CellDataVectorCreator
        Assert ((dynamic_cast<DataPostprocessor<dim>*>(viz_postprocessor)
                 != nullptr)
                ||
                (dynamic_cast<VisualizationPostprocessors::CellDataVectorCreator<dim>*>(viz_postprocessor)
                 != nullptr)
                ,
                ExcMessage ("Can't convert visualization postprocessor to type "
                            "dealii::DataPostprocessor or "
                            "VisualizationPostprocessors::CellDataVectorCreator!?"));

        postprocessors.push_back (std::unique_ptr<VisualizationPostprocessors::Interface<dim> >
                                  (viz_postprocessor));

        if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(&*postprocessors.back()))
          sim->initialize_simulator (this->get_simulator());

        postprocessors.back()->parse_parameters (prm);
        postprocessors.back()->initialize ();
      }
    }


    template <int dim>
    void
    Visualization<dim>::
    register_visualization_postprocessor (const std::string &name,
                                          const std::string &description,
                                          void (*declare_parameters_function) (ParameterHandler &),
                                          VisualizationPostprocessors::Interface<dim> *(*factory_function) ())
    {
      std::get<dim>(registered_visualization_plugins).register_plugin (name,
                                                                       description,
                                                                       declare_parameters_function,
                                                                       factory_function);
    }


    template <int dim>
    std::list<std::string>
    Visualization<dim>::required_other_postprocessors () const
    {
      std::list<std::string> requirements;

      // loop over all of the viz postprocessors and collect what
      // they want. don't worry about duplicates, the postprocessor
      // manager will filter them out
      for (const auto &p : postprocessors)
        {
          const std::list<std::string> this_requirements = p->required_other_postprocessors();
          requirements.insert (requirements.end(),
                               this_requirements.begin(), this_requirements.end());
        }

      return requirements;
    }


    template <int dim>
    void
    Visualization<dim>::set_last_output_time (const double current_time)
    {
      // if output_interval is positive, then update the last supposed output
      // time
      if (output_interval > 0)
      {
        // We need to find the last time output was supposed to be written.
        // this is the last_output_time plus the largest positive multiple
        // of output_intervals that passed since then. We need to handle the
        // edge case where last_output_time+output_interval==current_time,
        // we did an output and std::floor sadly rounds to zero. This is done
        // by forcing std::floor to round 1.0-eps to 1.0.
        const double magic = 1.0+2.0*std::numeric_limits<double>::epsilon();
        last_output_time = last_output_time + std::floor((current_time-last_output_time)/output_interval*magic) * output_interval/magic;
      }
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace internal
  {
    namespace Plugins
    {
      template <>
      std::list<internal::Plugins::PluginList<Postprocess::VisualizationPostprocessors::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<Postprocess::VisualizationPostprocessors::Interface<2> >::plugins = nullptr;
      template <>
      std::list<internal::Plugins::PluginList<Postprocess::VisualizationPostprocessors::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<Postprocess::VisualizationPostprocessors::Interface<3> >::plugins = nullptr;
    }
  }

  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
#define INSTANTIATE(dim) \
  template class Interface<dim>;

      ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }


  namespace Postprocess
  {
    ELASPECT_REGISTER_POSTPROCESSOR(Visualization,
                                "visualization",
                                "A postprocessor that takes the solution and writes "
                                "it into files that can be read by a graphical "
                                "visualization program. Additional run time parameters "
                                "are read from the parameter subsection 'Visualization'.")
  }
}
