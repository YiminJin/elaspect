#include <elaspect/structured_data.h>
#include <elaspect/utilities.h>
#include <elaspect/geometry_model/box.h>
#include <elaspect/geometry_model/sphere.h>
#include <elaspect/geometry_model/spherical_shell.h>

#include <deal.II/base/function_lib.h>

#include <boost/lexical_cast.hpp>

namespace elaspect
{
  namespace Utilities
  {
    /*--------------------- class StructuredDataLookup ------------------*/

    template <int dim>
    StructuredDataLookup<dim>::StructuredDataLookup(const unsigned int components,
                                                    const double scale_factor)
      :
      components(components),
      data_component_names(components),
      data(components),
      maximum_component_value(components),
      scale_factor(scale_factor),
      coordinate_values_are_equidistant(false)
    {}


    template <int dim>
    StructuredDataLookup<dim>::StructuredDataLookup(const double scale_factor)
      :
      components(numbers::invalid_unsigned_int),
      scale_factor(scale_factor),
      coordinate_values_are_equidistant(false)
    {}


    namespace
    {
      /**
       * Return whether a set of coordinate points in dim space dimensions
       * are all equidistantly spaced within some small tolerance.
       */
      template <int dim>
      bool data_is_equidistant(const std::array<std::vector<double>, dim> &coordinate_values)
      {
        bool coordinate_values_are_equidistant = true;

        for (unsigned int d = 0; d < dim; ++d)
        {
          const double grid_spacing = coordinate_values[d][1] - coordinate_values[d][0];

          for (unsigned int n = 1; n < coordinate_values[d].size(); ++n)
          {
            const double current_grid_spacing = coordinate_values[d][n] - coordinate_values[d][n-1];

            AssertThrow(current_grid_spacing > 0.0,
                        ExcMessage("Coordinates in dimension "
                                   + Utilities::int_to_string(d)
                                   + " are not strictly ascending."));

                // If spacing between coordinates changed (with a relative
                // tolerance), keep track of that information.  Note that we do
                // not break out of this loop in this case but run through the
                // whole array, so that the AssertThrow above is executed for
                // each entry to ensure increasing coordinate values.
                if (std::abs(current_grid_spacing - grid_spacing) > 
                    0.005 * (current_grid_spacing + grid_spacing))
                  coordinate_values_are_equidistant = false;
          }
        }

        return coordinate_values_are_equidistant;
      }
    }


    template <int dim>
    void
    StructuredDataLookup<dim>::reinit(const std::vector<std::string> &column_names,
                                      std::vector<std::vector<double>> &&coordinate_values_,
                                      std::vector<Table<dim,double>> &&data_table,
                                      const MPI_Comm &mpi_communicator,
                                      const unsigned int root_process)
    {
      // If this is the root process, or if the user did not request
      // sharing, then set up the various member variables we need
      // to compute from the input data
      if (root_process == numbers::invalid_unsigned_int ||
          Utilities::MPI::this_mpi_process(mpi_communicator) == root_process)
      {
        Assert(coordinate_values_.size() == dim,
               ExcMessage("Invalid size of coordinate_values."));
        for (unsigned int d = 0; d < dim; ++d)
        {
          this->coordinate_values[d] = std::move(coordinate_values_[d]);
          AssertThrow(this->coordinate_values[d].size() > 1,
                      ExcMessage("Error: At least 2 entries per coordinate direction are required"));
          table_points[d] = this->coordinate_values[d].size();
        }

        components = column_names.size();
        data_component_names = column_names;
        Assert(data_table.size() == components,
               ExcMessage("Error: Incorrect number of columns specified."));
        for (unsigned int c = 0; c < components; ++c)
          Assert(data_table[c].size() == table_points,
                 ExcMessage("Error: One of the data tables has an incorrect size."));

        // compute maximum_component_value for each component:
        maximum_component_value = std::vector<double>(components, 
                                                      -std::numeric_limits<double>::max());
        for (unsigned int c = 0; c < components; ++c)
        {
          const std::size_t n_elements = data_table[c].n_elements();
          for (std::size_t idx = 0; idx < n_elements; ++idx)
            maximum_component_value[c] = 
              std::max(maximum_component_value[c], 
                       data_table[c](compute_table_indices(table_points, idx)));
        }

        // In case the data is specified on a grid that is equidistant
        // in each coordinate direction, we only need to store
        // (besides the data) the number of intervals in each direction and
        // the begin- and end-points of the coordinates.
        // In case the grid is not equidistant, we need to keep
        // all the coordinates in each direction, which is more costly.
        coordinate_values_are_equidistant = data_is_equidistant<dim> (coordinate_values);
      }

      // If the caller of this function has actually requested sharing data,
      // then we have set up member variables on the root process, but not
      // on any of the other processes. Broadcast the data to the remaining
      // processes
      if (root_process != numbers::invalid_unsigned_int)
      {
        coordinate_values                 = Utilities::MPI::broadcast (mpi_communicator,
                                                                       coordinate_values,
                                                                       root_process);
        components                        = Utilities::MPI::broadcast (mpi_communicator,
                                                                       components,
                                                                       root_process);
        data_component_names              = Utilities::MPI::broadcast (mpi_communicator,
                                                                       data_component_names,
                                                                       root_process);
        maximum_component_value           = Utilities::MPI::broadcast (mpi_communicator,
                                                                       maximum_component_value,
                                                                       root_process);
        coordinate_values_are_equidistant = Utilities::MPI::broadcast (mpi_communicator,
                                                                       coordinate_values_are_equidistant,
                                                                       root_process);
        table_points                      = Utilities::MPI::broadcast (mpi_communicator,
                                                                       table_points,
                                                                       root_process);

        // We can then also prepare the data tables for sharing between
        // processes
        for (unsigned int c = 0; c < components; ++c)
          data_table[c].replicate_across_communicator (mpi_communicator,
                                                       root_process);
      }

      Assert(data_table.size() == components,
             ExcMessage("Error: Incorrect number of columns specified."));
      for (unsigned int c = 0; c < components; ++c)
        Assert(data_table[c].size() == table_points,
               ExcMessage("Error: One of the data tables has an incorrect size."));

      // For each data component, set up a GridData,
      // its type depending on the read-in grid.
      data.resize(components);
      for (unsigned int c = 0; c < components; ++c)
      {
        if (coordinate_values_are_equidistant)
        {
          std::array<unsigned int,dim> table_intervals;
          for (unsigned int d = 0; d < dim; ++d)
            table_intervals[d] = table_points[d]-1;

          // The min and max of the coordinates in the data file.
          std::array<std::pair<double,double>,dim> grid_extent;
          for (unsigned int d = 0; d < dim; ++d)
          {
            grid_extent[d].first = coordinate_values[d][0];
            grid_extent[d].second = coordinate_values[d][table_points[d]-1];

            Assert(table_intervals[d] >= 1,
                   ExcMessage("There needs to be at least one subinterval in each "
                              "coordinate direction."));
            Assert(grid_extent[d].first < grid_extent[d].second,
                   ExcMessage("The interval in each coordinate direction needs "
                              "to have positive size"));
          }

          data[c]
            = std::make_unique<Functions::InterpolatedUniformGridData<dim>>
              (std::move(grid_extent),
               std::move(table_intervals),
               std::move(data_table[c]));           
        }
        else
          // Create the object and move the big objects. Due to an old design flaw,
          // the current class stores a copy of the 'coordinate_values' and some
          // plugins actually use it too, i.e., we can't move the data out of
          // this object. In another design flaw, the deal.II classes
          // do not make the coordinate values accessible, so we have to continue
          // storing a copy. In other words, we can't just move stuff --
          // we have to make a copy of 'coordinate_values' and then move
          // that copy.
          //
          // (The call to std::move on the first argument is unnecessary: We
          // create a temporary object, and that's an rvalue that the constructor
          // we call would bind to. But never a bad idea to be explicit.)
          data[c]
            = std::make_unique<Functions::InterpolatedTensorProductGridData<dim>>
              (std::move(std::array<std::vector<double>,dim>(this->coordinate_values)),
               std::move(data_table[c]));
      }
    }


    template <int dim>
    void 
    StructuredDataLookup<dim>::load_file(const std::string &filename,
                                         const MPI_Comm &comm)
    {
      const unsigned int root_process = 0;

      std::vector<std::string> column_names;
      std::vector<Table<dim,double>> data_tables;
      std::vector<std::vector<double>> coordinate_values(dim);

      // If this is the root process, then set up the various member variables
      // we need to compute from the input data
      if (Utilities::MPI::this_mpi_process(comm) == root_process)
      {
        // Grab the values already stored in this class (if they exist), this way we can
        // check if somebody changes the size of the table over time and error out (see below)
        TableIndices<dim> new_table_points = this->table_points;

        // Read data from disk and distribute among processes. 
        std::stringstream in(read_and_distribute_file_content(filename, MPI_COMM_SELF));

        // Read header lines and table size
        while (in.peek() == '#')
        {
          std::string line;
          std::getline(in,line);
          std::stringstream linestream(line);
          std::string word;
          while (linestream >> word)
            if (word == "POINTS:")
              for (unsigned int i = 0; i < dim; i++)
              {
                unsigned int temp_index;
                linestream >> temp_index;

                if (new_table_points[i] == 0)
                  new_table_points[i] = temp_index;
                else
                  AssertThrow (new_table_points[i] == temp_index,
                               ExcMessage("The file grid must not change over model runtime. "
                                          "Either you prescribed a conflicting number of points in "
                                          "the input file, or the POINTS comment in your data files "
                                          "is changing between following files."));
              }
        }

        for (unsigned int i = 0; i < dim; i++)
        {
          AssertThrow(new_table_points[i] != 0,
                      ExcMessage("Could not successfully read in the file header of the "
                                 "ascii data file <" + filename + ">. One header line has to "
                                 "be of the format: '#POINTS: N1 [N2] [N3]', where N1 and "
                                 "potentially N2 and N3 have to be the number of data points "
                                 "in their respective dimension. Check for typos in this line "
                                 "(e.g. a missing space character)."));
        }

        // Read column lines if present
        unsigned int name_column_index = 0;
        double temp_data;

        while(true)
        {
          AssertThrow (name_column_index < 100,
                       ExcMessage("The program found more than 100 columns in the first line of the data file. "
                                  "This is unlikely intentional. Check your data file and make sure the data can be "
                                  "interpreted as floating point numbers. If you do want to read a data file with more "
                                  "than 100 columns, please remove this assertion."));

          std::string column_name_or_data;
          in >> column_name_or_data;
          try
          {
            // If the data field contains a name this will throw an exception
            temp_data = boost::lexical_cast<double>(column_name_or_data);

            // If there was no exception we have left the line containing names
            // and have read the first data field. Save number of components, and
            // make sure there is no contradiction if the components were already 
            // given to the constructor of this class.
            if (components == numbers::invalid_unsigned_int)
              components = name_column_index - dim;
            else if (name_column_index != 0)
              AssertThrow(components + dim == name_column_index,
                          ExcMessage("The number of expected data columns and the "
                                     "list of column names at the beginning of the data file "
                                     + filename + " do not match. The file should contain "
                                     + Utilities::int_to_string(name_column_index) + " column "
                                     "names (one for each dimension and one per data column), "
                                     "but it only has " + Utilities::int_to_string(components+dim) +
                                     " column names."));
            break;
          }
          catch (const boost::bad_lexical_cast &e)
          {
            // The first dim columns are coordinates and contain no data
            if (name_column_index >= dim)
            {
              // Transform name to lower case to prevent confusion with capital letters
              // Note: only ASCII characters allowed
              std::transform(column_name_or_data.begin(),
                             column_name_or_data.end(),
                             column_name_or_data.begin(),
                             ::tolower);

              AssertThrow(std::find(column_names.begin(), 
                                    column_names.end(), 
                                    column_name_or_data) 
                          == column_names.end(),
                          ExcMessage("There are multiple fields named " + column_name_or_data +
                                     " in the data file " + filename + ". Please remove duplication to "
                                     "allow for unique association between column and name."));

              column_names.push_back(column_name_or_data);
            }
            ++name_column_index;
          }
        }

        // Create table for the data. This peculiar reinit is necessary, because
        // there is no constructor for Table that takes TableIndices as
        // argument.
        Table<dim,double> data_table;
        data_table.TableBase<dim,double>::reinit(new_table_points);
        AssertThrow (components != numbers::invalid_unsigned_int,
                     ExcMessage("ERROR: number of components in " + filename + " could not be "
                                "determined automatically. Either add a header with column "
                                "names or pass the number of columns in the StructuredData "
                                "constructor."));
        data_tables.resize(components, data_table);

        for (unsigned int d = 0; d < dim; ++d)
          coordinate_values[d].resize(new_table_points[d]);

        if (column_names.size() == 0)
        {
          // set default column names:
          for (unsigned int c=0; c<components; ++c)
            column_names.push_back("column " + Utilities::int_to_string(c,2));
        }

        // Make sure that data file actually has as many columns as we think it has
        // (either based on the header, or based on what was passed to the constructor).
        const std::streampos position = in.tellg();
        std::string first_data_row;
        std::getline(in, first_data_row);
        std::stringstream linestream(first_data_row);
        std::string column_entry;

        // We have already read in the first data entry above in the try/catch block,
        // so there's one more column in the file than in the line we just read in.
        unsigned int number_of_entries = 1;
        while (linestream >> column_entry)
          number_of_entries += 1;

        AssertThrow ((number_of_entries) == column_names.size()+dim,
                     ExcMessage("ERROR: The number of columns in the data file " + filename +
                                " is incorrect. It needs to have " + Utilities::int_to_string(column_names.size()+dim) +
                                " columns, but the first row has " + Utilities::int_to_string(number_of_entries) +
                                " columns."));

        // Go back to the position in the file where we started the check for the column numbers.
        in.seekg (position);

        // Finally read data lines:
        std::size_t read_data_entries = 0;

        do
        {
          // what row and column of the file are we in?
          const std::size_t column_num = read_data_entries%(components+dim);
          const std::size_t row_num = read_data_entries/(components+dim);
          const TableIndices<dim> idx = compute_table_indices(new_table_points, row_num);

          if (column_num < dim)
          {
            // This is a coordinate. Store (and check that they are consistent)
            const double old_value = coordinate_values[column_num][idx[column_num]];

            AssertThrow(old_value == 0. ||
                        (std::abs(old_value-temp_data) < 1e-8*std::abs(old_value)),
                        ExcMessage("Invalid coordinate in column "
                                   + Utilities::int_to_string(column_num) + " in row "
                                   + Utilities::int_to_string(row_num)
                                   + " in file " + filename +
                                   "\nThis class expects the coordinates to be structured, meaning "
                                   "the coordinate values in each coordinate direction repeat exactly "
                                   "each time. This also means each row in the data file has to have "
                                   "the same number of columns as the first row containing data."));

            coordinate_values[column_num][idx[column_num]] = temp_data;
          }
          else
          {
            // This is a data value, so scale and store:
            const unsigned int component = column_num - dim;
            data_tables[component](idx) = temp_data * scale_factor;
          }

          ++read_data_entries;
        }
        while (in >> temp_data);

        AssertThrow(in.eof(),
                    ExcMessage ("While reading the data file '" + filename + "' the ascii data "
                                "plugin has encountered an error before the end of the file. "
                                "Please check for malformed data values (e.g. NaN) or superfluous "
                                "lines at the end of the data file."));

        const std::size_t n_expected_data_entries = (components + dim) * data_table.n_elements();
        AssertThrow(read_data_entries == n_expected_data_entries,
                    ExcMessage ("While reading the data file '" + filename + "' the ascii data "
                                "plugin has reached the end of the file, but has not found the "
                                "expected number of data values considering the spatial dimension, "
                                "data columns, and number of lines prescribed by the POINTS header "
                                "of the file. Please check the number of data "
                                "lines against the POINTS header in the file."));
      }

      // We have set up member variables on the root process, but not on any of
      // the other processes. So broadcast the data to the remaining processes.
      components = Utilities::MPI::broadcast (comm,
                                              components,
                                              root_process);
      coordinate_values = Utilities::MPI::broadcast (comm,
                                                     coordinate_values,
                                                     root_process);
      column_names = Utilities::MPI::broadcast (comm,
                                                column_names,
                                                root_process);

      if (Utilities::MPI::this_mpi_process(comm) != root_process)
        data_tables.resize(components);

      // Finally create the data. We want to call the move-version of reinit() so
      // that the data doesn't have to be copied, so use std::move on all big
      // objects.
      this->reinit(column_names,
                   std::move(coordinate_values),
                   std::move(data_tables),
                   comm,
                   root_process);
    }


    template <int dim>
    TableIndices<dim>
    StructuredDataLookup<dim>::compute_table_indices(const TableIndices<dim> &sizes,
                                                     const std::size_t idx) const
    {
      TableIndices<dim> result;
      result[0] = idx % sizes[0];
      if (dim >= 2)
        result[1] = (idx / sizes[0]) % sizes[1];
      if (dim == 3)
        result[2] = idx / (sizes[0] * sizes[1]);

      return result;
    }


    template <int dim>
    double
    StructuredDataLookup<dim>::get_data(const Point<dim> &position,
                                        const unsigned int component) const
    {
      Assert(component<components, ExcMessage("Invalid component index"));
      return data[component]->value(position);
    }


    template <int dim>
    Tensor<1,dim>
    StructuredDataLookup<dim>::get_gradients(const Point<dim> &position,
                                             const unsigned int component) const
    {
      return data[component]->gradient(position,0);
    }


    template <int dim>
    double
    StructuredDataLookup<dim>::get_maximum_component_value(const unsigned int component) const
    {
      return maximum_component_value[component];
    }


    /*------------------------ class AsciiDataBase ----------------------*/

    template <int dim>
    AsciiDataBase<dim>::AsciiDataBase() 
     = default;


    template <int dim>
    void
    AsciiDataBase<dim>::declare_parameters (ParameterHandler  &prm,
                                            const std::string &default_directory,
                                            const std::string &default_filename,
                                            const std::string &subsection_name)
    {
      prm.enter_subsection (subsection_name);
      {
        prm.declare_entry ("Data directory",
                           default_directory,
                           Patterns::DirectoryName (),
                           "The name of a directory that contains the model data. This path "
                           "may either be absolute (if starting with a `/') or relative to "
                           "the current directory. The path may also include the special "
                           "text `$ELASPECT_SOURCE_DIR' which will be interpreted as the path "
                           "in which the elASPECT source files were located when elASPECT was "
                           "compiled. This interpretation allows, for example, to reference "
                           "files located in the `data/' subdirectory of elASPECT.");
        prm.declare_entry ("Data file name",
                           default_filename,
                           Patterns::Anything (),
                           "The file name of the model data.");
        prm.declare_entry ("Scale factor", "1.",
                           Patterns::Double (),
                           "Scalar factor, which is applied to the model data. "
                           "You might want to use this to scale the input to a "
                           "reference model. Another way to use this factor is to "
                           "convert units of the input files. For instance, if you "
                           "provide velocities in cm/yr set this factor to 0.01.");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    AsciiDataBase<dim>::parse_parameters (ParameterHandler &prm,
                                          const std::string &subsection_name)
    {
      prm.enter_subsection (subsection_name);
      {
        // Get the path to the data files. If it contains a reference
        // to $ELASPECT_SOURCE_DIR, replace it by what CMake has given us
        // as a #define
        data_directory = Utilities::expand_ELASPECT_SOURCE_DIR(prm.get ("Data directory"));
        data_file_name    = prm.get ("Data file name");
        scale_factor      = prm.get_double ("Scale factor");
      }
      prm.leave_subsection();
    }


    /*---------------------- class AsciiDataBoundary --------------------*/

    namespace
    {
      template <int dim>
      void
      check_supported_geometry_models(const GeometryModel::Interface<dim> &geometry_model)
      {
        AssertThrow ((Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>>(geometry_model))
                     || (Plugins::plugin_type_matches<const GeometryModel::Sphere<dim>>(geometry_model))
                     || (Plugins::plugin_type_matches<const GeometryModel::Box<dim>>(geometry_model)),
                     ExcMessage ("This ascii data plugin can only be used with supported "
                                 "geometry models."));
      }
    }


    template <int dim>
    AsciiDataBoundary<dim>::AsciiDataBoundary()
      :
      current_file_number(numbers::invalid_unsigned_int),
      first_data_file_number(numbers::invalid_unsigned_int),
      time_weight(numbers::signaling_nan<double>()),
      time_dependent(false),
      lookups(),
      old_lookups()
    {}


    template <int dim>
    void
    AsciiDataBoundary<dim>::initialize(const std::set<types::boundary_id> &boundary_ids,
                                       const unsigned int components)
    {
      check_supported_geometry_models(this->get_geometry_model());

      for (const auto &boundary_id : boundary_ids)
      {
        lookups.insert(std::make_pair(boundary_id,
                                      std::make_unique<Utilities::StructuredDataLookup<dim-1>>
                                      (components,
                                       this->scale_factor)));

        // Set the first file number and load the first files
        current_file_number = first_data_file_number;

        const std::string filename(create_filename(current_file_number, boundary_id));

        this->get_pcout() << std::endl << "   Loading Ascii data boundary file "
                          << filename << '.' << std::endl << std::endl;

        AssertThrow(Utilities::fexists(filename, this->get_mpi_communicator()),
                    ExcMessage(std::string("Ascii data file <") + filename + "> not found!"));

        lookups.find(boundary_id)->second->load_file(filename, this->get_mpi_communicator());

        if (time_dependent)
        {
          old_lookups.insert(std::make_pair(boundary_id,
                                            std::make_unique<Utilities::StructuredDataLookup<dim-1>>
                                            (components,
                                             this->scale_factor)));

          const std::string filename (create_filename (current_file_number + 1, boundary_id));
          if (Utilities::fexists(filename, this->get_mpi_communicator()))
          {
            this->get_pcout() << std::endl << "   Also loading next Ascii data boundary file "
                              << filename << '.' << std::endl << std::endl;
            lookups.find(boundary_id)->second.swap(old_lookups.find(boundary_id)->second);
            lookups.find(boundary_id)->second->load_file(filename, this->get_mpi_communicator());
          }
          else
          {
            // next file not found, end looking for new files
            time_dependent = false;
          }
        }
      }
    }


    namespace
    {
      /**
       * Given a string @p filename_and_path that contains exactly one
       * <code>%s</code> and one <code>%d</code> code (possibly modified
       * by flag, field, and length modifiers as discussed in the man
       * pages of the <code>printf()</code> family of functions),
       * return the expanded string where the <code>%s</code> code is
       * replaced by @p boundary_name, and <code>%d</code> is replaced
       * by @p filenumber.
       */
      std::string replace_placeholders(const std::string &filename_and_path,
                                       const std::string &boundary_name,
                                       const int filenumber)
      {
        const int maxsize = filename_and_path.length() + 256;
        char *filename = static_cast<char *>(malloc (maxsize * sizeof(char)));
        int ret = snprintf (filename,
                            maxsize,
                            filename_and_path.c_str(),
                            boundary_name.c_str(),
                            filenumber);

        AssertThrow(ret >= 0, ExcMessage("Invalid string placeholder in filename detected."));
        AssertThrow(ret< maxsize, ExcInternalError("snprintf string overflow detected."));
        const std::string str_result (filename);
        free (filename);
        return str_result;
      }
    }


    template <int dim>
    std::string
    AsciiDataBoundary<dim>::create_filename(const int filenumber,
                                            const types::boundary_id boundary_id) const
    {
      std::string templ = this->data_directory + this->data_file_name;

      const std::string boundary_name = this->get_geometry_model().translate_id_to_symbol_name(boundary_id);

      const std::string result = replace_placeholders(templ, boundary_name, filenumber);
      if (fexists(result, this->get_mpi_communicator()))
        return result;

      // Backwards compatibility check: people might still be using the old
      // names of the top/bottom boundary. If they do, print a warning but
      // accept those files.
      std::string compatible_result;
      if (boundary_name == "top")
      {
        compatible_result = replace_placeholders(templ, "surface", filenumber);
        if (!fexists(compatible_result, this->get_mpi_communicator()))
          compatible_result = replace_placeholders(templ, "outer", filenumber);
      }
      else if (boundary_name == "bottom")
        compatible_result = replace_placeholders(templ, "inner", filenumber);

      if (!fexists(result, this->get_mpi_communicator()) && fexists(compatible_result, this->get_mpi_communicator()))
      {
        this->get_pcout() << "WARNING: Filename convention concerning geometry boundary "
                          "names changed. Please rename '" << compatible_result << "'"
                          << " to '" << result << "'"
                          << std::endl;
        return compatible_result;
      }

      return result;
    }


    namespace
    {
      /**
       * Determines which of the dimensions of a position are aligned with
       * a certain @p boundary_id. E.g. the left boundary of a box
       * model extents in the y and z direction (position[1] and
       * position[2]), therefore the function would return [1,2] for dim==3
       * or [1] for dim==2. We are lucky that these indices are identical
       * for all existing geometries (if we use natural coordinates),
       * therefore we do not need to distinguish between them. If we
       * introduce a geometry model for which boundaries are not aligned
       * with natural coordinates this needs to become a function in
       * the interface of the geometry model and needs to be specialized
       * for distinct geometry models.
       */
      template <int dim>
      std::array<unsigned int, dim-1>
      get_boundary_dimensions(const types::boundary_id boundary_id)
      {
        std::array<unsigned int,dim-1> boundary_dimensions;
        boundary_dimensions.fill(numbers::invalid_unsigned_int);

        switch (dim)
        {
          case 2:
            if ((boundary_id == 2) || (boundary_id == 3) || (boundary_id == 4) || (boundary_id == 5))
            {
              boundary_dimensions[0] = 0;
            }
          else if ((boundary_id == 0) || (boundary_id == 1))
            {
              boundary_dimensions[0] = 1;
            }
            else
            {
              AssertThrow(false,ExcNotImplemented());
            }
            break;

          case 3:
            if ((boundary_id == 4) || (boundary_id == 5) || (boundary_id == 8) || (boundary_id == 9))
            {
              boundary_dimensions[0] = 0;
              boundary_dimensions[1] = 1;
            }
            else if ((boundary_id == 0) || (boundary_id == 1))
            {
              boundary_dimensions[0] = 1;
              boundary_dimensions[1] = 2;
            }
            else if ((boundary_id == 2) || (boundary_id == 3) || (boundary_id == 6) || (boundary_id == 7))
            {
              boundary_dimensions[0] = 0;
              boundary_dimensions[1] = 2;
            }
            else
            {
              AssertThrow(false,ExcNotImplemented());
            }
            break;

          default:
            AssertThrow(false,ExcNotImplemented());
        }

        return boundary_dimensions;
      }


      /**
       * Convert a point @p position in cartesian coordinates into the set of coordinates
       * used in the structured data input files, i.e. coordinates in the
       * natural coordinate system of the given @p geometry_model.
       */
      template <int dim>
      Point<dim>
      data_coordinates_from_position(const Point<dim> &position,
                                     const GeometryModel::Interface<dim> &geometry_model)
      {
        const std::array<double,dim> natural_position = geometry_model.cartesian_to_natural_coordinates(position);
        Point<dim> data_coordinates = Utilities::convert_array_to_point<dim>(natural_position);

        return data_coordinates;
      }


      /**
       * Convert a point @p data_coordinates using coordinates appropriate for
       * volumetric structured data input files (e.g. produced by the
       * data_coordinates_from_position() function above) into appropriate
       * coordinates on the boundary given by @p boundary_indicator.
       * Because all existing geometry models use boundaries that are aligned
       * with natural coordinate directions this means finding the correct components
       * of the input coordinates and returning a point with one dimension less
       * than the input point that is determined by these coordinates.
       */
      template <int dim>
      Point<dim-1>
      boundary_coordinates_from_data_coordinates(const Point<dim> &data_coordinates,
                                                 const types::boundary_id boundary_indicator)
      {
        const std::array<unsigned int,dim-1> boundary_dimensions =
          get_boundary_dimensions<dim>(boundary_indicator);

        Point<dim-1> boundary_coordinates;
        for (unsigned int i = 0; i < dim-1; ++i)
          boundary_coordinates[i] = data_coordinates[boundary_dimensions[i]];

        return boundary_coordinates;
      }
    }


    template <int dim>
    double
    AsciiDataBoundary<dim>::
    get_data_component(const types::boundary_id boundary_indicator,
                       const Point<dim> &       position,
                       const unsigned int       component) const
    {
      const Point<dim> data_coordinates = data_coordinates_from_position(position, this->get_geometry_model());
      const Point<dim-1> boundary_coordinates = boundary_coordinates_from_data_coordinates(data_coordinates, boundary_indicator);

      Assert(lookups.find(boundary_indicator) != lookups.end(),
             ExcInternalError());
      const double data = lookups.find(boundary_indicator)->second->get_data(boundary_coordinates, component);

      if (!time_dependent)
        return data;

      const double old_data = old_lookups.find(boundary_indicator)->second->get_data(boundary_coordinates, component);

      return time_weight * data + (1 - time_weight) * old_data;
    }


    template <int dim>
    double
    AsciiDataBoundary<dim>::
    get_maximum_component_value(const types::boundary_id boundary_indicator,
                                const unsigned int component) const
    {
      return lookups.find(boundary_indicator)->second->get_maximum_component_value(component);
    }


    // explicit instantiations
    template class StructuredDataLookup<1>;
    template class StructuredDataLookup<2>;
    template class StructuredDataLookup<3>;
    template class AsciiDataBase<2>;
    template class AsciiDataBase<3>;
    template class AsciiDataBoundary<2>;
    template class AsciiDataBoundary<3>;
  }
}
