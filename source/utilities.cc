#include <elaspect/utilities.h>
#include <elaspect/geometry_model/interface.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/lac/vector.h>

#include <dirent.h>
#include <fstream>
#include <sys/stat.h>
#include <regex>

#include <boost/iostreams/copy.hpp>

namespace elaspect
{
  namespace Utilities
  {
    bool fexists(const std::string &filename)
    {
      std::ifstream ifile(filename.c_str());

      // return whether construction of the input file has succeeded;
      // success requires the file to exist and to be readable
      return static_cast<bool>(ifile);
    }


    bool fexists(const std::string &filename,
                 MPI_Comm comm)
    {
      bool file_exists = false;
      if (Utilities::MPI::this_mpi_process(comm) == 0)
      {
        std::ifstream ifile(filename.c_str());

        // return whether construction of the input file has succeeded;
        // success requires the file to exist and to be readable
        file_exists = static_cast<bool>(ifile);
      }
      return Utilities::MPI::broadcast(comm, file_exists);
    }


    std::string
    read_and_distribute_file_content(const std::string &filename,
                                     const MPI_Comm &comm)
    {
      std::string data_string;

      if (Utilities::MPI::this_mpi_process(comm) == 0)
      {
        std::size_t filesize;
        std::ifstream filestream;

        const bool filename_ends_in_gz = std::regex_search(filename, std::regex("\\.gz$"));
        if (filename_ends_in_gz)
          filestream.open(filename.c_str(), std::ios_base::in | std::ios_base::binary);
        else
          filestream.open(filename.c_str());

        if (!filestream)
        {
          // broadcast failure state, then throw
          std::size_t invalid_filesize = numbers::invalid_size_type;
          const int ierr = MPI_Bcast(&invalid_filesize, 
                                     1,
                                     Utilities::MPI::mpi_type_id_for_type<std::size_t>, 
                                     0, 
                                     comm);
          AssertThrowMPI(ierr);
          AssertThrow(false,
                      ExcMessage(std::string("Could not open file <") + filename + ">."));
        }

        // Read data from disk
        std::stringstream datastream;

        try
        {
          boost::iostreams::filtering_istreambuf in;
          if (filename_ends_in_gz)
            in.push(boost::iostreams::gzip_decompressor());

          in.push(filestream);
          boost::iostreams::copy(in, datastream);
        }
        catch (const std::ios::failure &)
        {
          // broadcast failure state, then throw
          std::size_t invalid_filesize = numbers::invalid_size_type;
          const int ierr = MPI_Bcast(&invalid_filesize,
                                     1,
                                     Utilities::MPI::mpi_type_id_for_type<std::size_t>,
                                     0,
                                     comm);
          AssertThrowMPI(ierr);
          AssertThrow(false,
                      ExcMessage(std::string("Could not read file content from <") + filename + ">."));
        }

        data_string = datastream.str();
        filesize = data_string.size();

        // Distribute data_size and data across processes
        const int ierr = MPI_Bcast(&filesize,
                                   1,
                                   Utilities::MPI::mpi_type_id_for_type<std::size_t>,
                                   0,
                                   comm);
        AssertThrowMPI(ierr);

        dealii::Utilities::MPI::broadcast(&data_string[0], filesize, 0, comm);
      }
      else
      {
        // Prepare for receiving data
        std::size_t filesize;
        const int ierr = MPI_Bcast(&filesize, 
                                   1, 
                                   Utilities::MPI::mpi_type_id_for_type<std::size_t>,
                                   0, 
                                   comm);
        AssertThrowMPI(ierr);
        if (filesize == numbers::invalid_size_type)
          throw QuietException();

        data_string.resize(filesize);

        // Receive and store data
        dealii::Utilities::MPI::broadcast(&data_string[0], filesize, 0, comm);
      }

      return data_string;
    }


    void
    collect_and_write_file_content(const std::string &filename,
                                   const std::string &file_content,
                                   const MPI_Comm &comm)
    {
      const std::vector<std::string> collected_content = Utilities::MPI::gather(comm, file_content);

      if (Utilities::MPI::this_mpi_process(comm) == 0)
      {
        std::ofstream filestream;
        filestream.open(filename.c_str());

        AssertThrow(filestream.good(),
                    ExcMessage(std::string("Could not open file <") + filename + ">."));

        try
        {
          for (const auto &content : collected_content)
            filestream << content;

          bool success = filestream.good();
          const int ierr = MPI_Bcast(&success, 1, Utilities::MPI::mpi_type_id_for_type<bool>, 0, comm);
          AssertThrowMPI(ierr);
        }
        catch (const std::ios::failure &)
        {
          // broadcast failure state, then throw
          bool success = false;
          const int ierr = MPI_Bcast(&success, 1, Utilities::MPI::mpi_type_id_for_type<bool>, 0, comm);
          AssertThrowMPI(ierr);
          AssertThrow(false,
                      ExcMessage(std::string("Could not write content to file <") + filename + ">."));
        }

        filestream.close();
      }
      else
      {
        // Check if the file has been written successfully
        bool success;
        int ierr = MPI_Bcast(&success, 1, Utilities::MPI::mpi_type_id_for_type<bool>, 0, comm);
        AssertThrowMPI(ierr);
        if (!success)
          throw QuietException();
      }
    }


    int
    mkdirp(std::string pathname,const mode_t mode)
    {
      // force trailing / so we can handle everything in loop
      if (pathname[pathname.size()-1] != '/')
        {
          pathname += '/';
        }

      size_t pre = 0;
      size_t pos;

      while ((pos = pathname.find_first_of('/',pre)) != std::string::npos)
        {
          const std::string subdir = pathname.substr(0,pos++);
          pre = pos;

          // if leading '/', first string is 0 length
          if (subdir.size() == 0)
            continue;

          int mkdir_return_value;
          if ((mkdir_return_value = mkdir(subdir.c_str(),mode)) && (errno != EEXIST))
            return mkdir_return_value;

        }

      return 0;
    }

    void create_directory(const std::string &pathname,
                          const MPI_Comm &comm,
                          bool silent)
    {
      // verify that the output directory actually exists. if it doesn't, create
      // it on processor zero
      int error;

      if ((Utilities::MPI::this_mpi_process(comm) == 0))
        {
          DIR *output_directory = opendir(pathname.c_str());
          if (output_directory == nullptr)
            {
              if (!silent)
                std::cout << "\n"
                          << "-----------------------------------------------------------------------------\n"
                          << "The output directory <" << pathname
                          << "> provided in the input file appears not to exist.\n"
                          << "elASPECT will create it for you.\n"
                          << "-----------------------------------------------------------------------------\n\n"
                          << std::endl;

              error = Utilities::mkdirp(pathname, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);

            }
          else
            {
              error = closedir(output_directory);
            }
          // Broadcast error code
          const int ierr = MPI_Bcast (&error, 1, MPI_INT, 0, comm);
          AssertThrowMPI(ierr);
          AssertThrow (error == 0,
                       ExcMessage (std::string("Can't create the output directory at <") + pathname + ">"));
        }
      else
        {
          // Wait to receive error code, and throw QuietException if directory
          // creation has failed
          const int ierr = MPI_Bcast (&error, 1, MPI_INT, 0, comm);
          AssertThrowMPI(ierr);
          if (error!=0)
            throw elaspect::QuietException();
        }
    }


    std::string
    expand_ELASPECT_SOURCE_DIR (const std::string &location)
    {
      return Utilities::replace_in_string(location,
                                          "$ELASPECT_SOURCE_DIR",
                                          ELASPECT_SOURCE_DIR);
    }


    bool
    has_unique_entries (const std::vector<std::string> &strings)
    {
      const std::set<std::string> set_of_strings(strings.begin(),strings.end());
      return (set_of_strings.size() == strings.size());
    }


    namespace Coordinates
    {
      template <int dim>
      std::array<double,dim>
      cartesian_to_spherical_coordinates(const Point<dim> &position)
      {
        std::array<double,dim> scoord;

        scoord[0] = position.norm(); // R
        scoord[1] = std::atan2(position(1),position(0)); // Phi
        if (scoord[1] < 0.0)
          scoord[1] += 2.0*numbers::PI; // correct phi to [0,2*pi]
        if (dim==3)
          {
            if (scoord[0] > std::numeric_limits<double>::min())
              scoord[2] = std::acos(position(2)/scoord[0]);
            else
              scoord[2] = 0.0;
          }
        return scoord;
      }


      template <int dim>
      Point<dim>
      spherical_to_cartesian_coordinates(const std::array<double,dim> &scoord)
      {
        Point<dim> ccoord;

        switch (dim)
          {
            case 2:
            {
              ccoord[0] = scoord[0] * std::cos(scoord[1]); // X
              ccoord[1] = scoord[0] * std::sin(scoord[1]); // Y
              break;
            }
            case 3:
            {
              ccoord[0] = scoord[0] * std::sin(scoord[2]) * std::cos(scoord[1]); // X
              ccoord[1] = scoord[0] * std::sin(scoord[2]) * std::sin(scoord[1]); // Y
              ccoord[2] = scoord[0] * std::cos(scoord[2]); // Z
              break;
            }
            default:
              Assert (false, ExcNotImplemented());
              break;
          }

        return ccoord;
      }


      template <int dim>
      Tensor<1, dim>
      spherical_to_cartesian_vector(const Tensor<1, dim> &spherical_vector,
                                    const Point<dim> &position)
      {
        Tensor<1, dim> cartesian_vector;

        const std::array<double, dim> r_phi_theta = cartesian_to_spherical_coordinates(position);

        switch (dim)
          {
            case 2:
            {
              const double phi = r_phi_theta[1];

              const double u_r   = spherical_vector[0];
              const double u_phi = spherical_vector[1];

              cartesian_vector[0] = std::cos(phi)*u_r
                                    - std::sin(phi)*u_phi; // X
              cartesian_vector[1] = std::sin(phi)*u_r
                                    + std::cos(phi)*u_phi; // Y

              break;
            }
            case 3:
            {
              const double phi   = r_phi_theta[1];
              const double theta = r_phi_theta[2];

              const double u_r     = spherical_vector[0];
              const double u_phi   = spherical_vector[1];
              const double u_theta = spherical_vector[2];

              cartesian_vector[0] = std::cos(phi)*std::sin(theta)*u_r
                                    - std::sin(phi)*u_phi
                                    + std::cos(phi)*std::cos(theta)*u_theta; // X
              cartesian_vector[1] = std::sin(phi)*std::sin(theta)*u_r
                                    + std::cos(phi)*u_phi
                                    + std::sin(phi)*std::cos(theta)*u_theta; // Y
              cartesian_vector[2] = std::cos(theta)*u_r
                                    - std::sin(theta)*u_theta; // Z
              break;
            }

            default:
              Assert (false, ExcNotImplemented());
              break;
          }

        return cartesian_vector;
      }


      CoordinateSystem
      string_to_coordinate_system(const std::string &coordinate_system)
      {
        if (coordinate_system == "cartesian")
          return cartesian;
        else if (coordinate_system == "spherical")
          return spherical;
        else if (coordinate_system == "depth")
          return Coordinates::depth;
        else
          AssertThrow(false, ExcNotImplemented());

        return Coordinates::invalid;
      }
    }


    template <int dim>
    Point<dim> convert_array_to_point(const std::array<double,dim> &array)
    {
      Point<dim> point;
      for (unsigned int i = 0; i < dim; i++)
        point[i] = array[i];

      return point;
    }


    template <int dim>
    std::array<double,dim> convert_point_to_array(const Point<dim> &point)
    {
      std::array<double,dim> array;
      for (unsigned int i = 0; i < dim; i++)
        array[i] = point[i];

      return array;
    }


    template <int dim>
    NaturalCoordinate<dim>::NaturalCoordinate(Point<dim> &position,
                                              const GeometryModel::Interface<dim> &geometry_model)
    {
      coordinate_system = geometry_model.natural_coordinate_system();
      coordinates = geometry_model.cartesian_to_natural_coordinates(position);
    }


    template <int dim>
    NaturalCoordinate<dim>::NaturalCoordinate(const std::array<double, dim> &coord,
                                              const Utilities::Coordinates::CoordinateSystem &coord_system) :
      coordinate_system (coord_system), coordinates (coord)
    {}


    template <int dim>
    const std::array<double,dim> &
    NaturalCoordinate<dim>::get_coordinates() const
    {
      return coordinates;
    }


    Operator::Operator()
      :
      op(uninitialized)
    {}


    Operator::Operator(const operation _op)
      :
      op(_op)
    {}


    double
    Operator::operator() (const double x, const double y) const
    {
      switch (op)
        {
          case Utilities::Operator::add:
          {
            return x + y;
          }
          case Utilities::Operator::subtract:
          {
            return x - y;
          }
          case Utilities::Operator::minimum:
          {
            return std::min(x,y);
          }
          case Utilities::Operator::maximum:
          {
            return std::max(x,y);
          }
          case Utilities::Operator::replace_if_valid:
          {
            if (std::isnan(y))
              return x;
            else
              return y;
          }
          default:
          {
            Assert (false, ExcInternalError());
          }
        }
      return numbers::signaling_nan<double>();
    }


    bool
    Operator::operator== (const operation other_op) const
    {
      return other_op == op;
    }


    std::vector<Operator> create_model_operator_list(const std::vector<std::string> &operator_names)
    {
      std::vector<Operator> operator_list(operator_names.size());
      for (unsigned int i=0; i<operator_names.size(); ++i)
        {
          // create operator list
          if (operator_names[i] == "add")
            operator_list[i] = Operator(Operator::add);
          else if (operator_names[i] == "subtract")
            operator_list[i] = Operator(Operator::subtract);
          else if (operator_names[i] == "minimum")
            operator_list[i] = Operator(Operator::minimum);
          else if (operator_names[i] == "maximum")
            operator_list[i] = Operator(Operator::maximum);
          else if (operator_names[i] == "replace if valid")
            operator_list[i] = Operator(Operator::replace_if_valid);
          else
            AssertThrow(false,
                        ExcMessage ("elASPECT only accepts the following operators: "
                                    "add, subtract, minimum, maximum, and replace if valid. But your parameter file "
                                    "contains: " + operator_names[i] + ". Please check your parameter file.") );
        }

      return operator_list;
    }


    const std::string get_model_operator_options()
    {
      return "add|subtract|minimum|maximum|replace if valid";
    }


    template <int dim>
    VectorFunctionFromVelocityFunctionObject<dim>::
    VectorFunctionFromVelocityFunctionObject
    (const unsigned int n_components,
     const std::function<Tensor<1,dim> (const Point<dim> &)> &function_object)
      :
      Function<dim>(n_components),
      function_object (function_object)
    {}



    template <int dim>
    double
    VectorFunctionFromVelocityFunctionObject<dim>::
    value (const Point<dim> &p,
           const unsigned int component) const
    {
      Assert (component < this->n_components,
              ExcIndexRange (component, 0, this->n_components));

      if (component < dim)
        {
          const Tensor<1,dim> v = function_object(p);
          return v[component];
        }
      else
        return 0;
    }



    template <int dim>
    void
    VectorFunctionFromVelocityFunctionObject<dim>::
    vector_value (const Point<dim>   &p,
                  Vector<double>     &values) const
    {
      AssertDimension(values.size(), this->n_components);

      // set everything to zero, and then the right components to their correct values
      values = 0;

      const Tensor<1,dim> v = function_object(p);
      for (unsigned int d=0; d<dim; ++d)
        values(d) = v[d];
    }
  }
}


// explicit instantiations
namespace elaspect
{
  namespace Utilities
  {
#define INSTANTIATE(dim) \
    template \
    Point<dim> Coordinates::spherical_to_cartesian_coordinates<dim>(const std::array<double,dim> &scoord); \
    \
    template \
    Tensor<1,dim> Coordinates::spherical_to_cartesian_vector<dim>(const Tensor<1,dim> &spherical_vector, \
                                                                  const Point<dim> &position); \
    \
    template \
      std::array<double,dim> Coordinates::cartesian_to_spherical_coordinates<dim>(const Point<dim> &position); \
    \
    template \
    Point<dim> convert_array_to_point<dim>(const std::array<double,dim> &array); \
    \
    template \
    std::array<double,dim> convert_point_to_array<dim>(const Point<dim> &point); \
    \
    template \
    class NaturalCoordinate<dim>; \
    \
    template \
    class VectorFunctionFromVelocityFunctionObject<dim>; \

    ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
