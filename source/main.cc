#include <elaspect/simulator.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/multithread_info.h>

#include <csignal>
#include <chrono>

#ifdef DEBUG
#ifdef ELASPECT_USE_FP_EXCEPTIONS
#include <cfenv>
#endif
#endif

#if ELASPECT_USE_SHARED_LIBS==1
#  include <dlfcn.h>
#  ifdef ELASPECT_HAVE_LINK_H
#    include <link.h>
#  endif
#endif


// get the value of a particular parameter from the contents of the input
// file. return an empty string if not found
std::string
get_last_value_of_parameter(const std::string &parameters,
                            const std::string &parameter_name)
{
  std::string return_value;

  std::istringstream x_file(parameters);
  while (x_file)
    {
      // get one line and strip spaces at the front and back
      std::string line;
      std::getline(x_file, line);
      while ((line.size() > 0) && (line[0] == ' ' || line[0] == '\t'))
        line.erase(0, 1);
      while ((line.size() > 0)
             && (line[line.size() - 1] == ' ' || line[line.size() - 1] == '\t'))
        line.erase(line.size() - 1, std::string::npos);
      // now see whether the line starts with 'set' followed by multiple spaces
      // if not, try next line
      if (line.size() < 4)
        continue;

      if ((line[0] != 's') || (line[1] != 'e') || (line[2] != 't')
          || !(line[3] == ' ' || line[3] == '\t'))
        continue;

      // delete the "set " and then delete more spaces if present
      line.erase(0, 4);
      while ((line.size() > 0) && (line[0] == ' ' || line[0] == '\t'))
        line.erase(0, 1);
      // now see whether the next word is the word we look for
      if (line.find(parameter_name) != 0)
        continue;

      line.erase(0, parameter_name.size());
      while ((line.size() > 0) && (line[0] == ' ' || line[0] == '\t'))
        line.erase(0, 1);

      // we'd expect an equals size here
      if ((line.size() < 1) || (line[0] != '='))
        continue;

      // remove comment
      std::string::size_type pos = line.find('#');
      if (pos != std::string::npos)
        line.erase (pos);

      // trim the equals sign at the beginning and possibly following spaces
      // as well as spaces at the end
      line.erase(0, 1);
      while ((line.size() > 0) && (line[0] == ' ' || line[0] == '\t'))
        line.erase(0, 1);
      while ((line.size() > 0) && (line[line.size()-1] == ' ' || line[line.size()-1] == '\t'))
        line.erase(line.size()-1, std::string::npos);

      // the rest should now be what we were looking for
      return_value = line;
    }

  return return_value;
}


unsigned int
get_dimension (const std::string &parameters)
{
  const std::string dimension = get_last_value_of_parameter(parameters, "Dimension");
  if (dimension.size() > 0)
  {
    AssertThrow (dimension.find('\r') == std::string::npos,
                 dealii::ExcMessage ("It appears that your input file uses Windows-style "
                                     "line endings ('\\r\\n') but you are running on a system where "
                                     "the C++ run time environment expects input files to have "
                                     "Unix-style line endings ('\\n'). You need to convert your "
                                     "input file to use the correct line endings before running "
                                     "elASPECT with it."));
    try
    {
      return dealii::Utilities::string_to_int (dimension);
    }
    catch (...)
    {
      AssertThrow (false,
                   dealii::ExcMessage("While reading the dimension from the input file, "
                                      "elASPECT found a string that can not be converted to "
                                      "an integer: <" + dimension + ">."));
      return 0; // we should never get here.
    }
  }
  else
    return 2;
}


#if ELASPECT_USE_SHARED_LIBS==1

#ifdef ELASPECT_HAVE_LINK_H
// collect the names of the shared libraries linked to by this program. this
// function is a callback for the dl_iterate_phdr() function we call below
int get_names_of_shared_libs (struct dl_phdr_info *info,
                              size_t,
                              void *data)
{
  reinterpret_cast<std::set<std::string>*>(data)->insert (info->dlpi_name);
  return 0;
}
#endif


// make sure the list of shared libraries we currently link with
// has deal.II only once
void validate_shared_lib_list (const bool before_loading_shared_libs)
{
#ifdef ELASPECT_HAVE_LINK_H
  // get the list of all shared libs we currently link against
  std::set<std::string> shared_lib_names;
  dl_iterate_phdr(get_names_of_shared_libs, &shared_lib_names);

  // find everything that is interesting
  std::set<std::string> dealii_shared_lib_names;
  for (const auto &p : shared_lib_names)
    if (p.find ("libdeal_II") != std::string::npos ||
        p.find ("libdeal.ii") != std::string::npos)
      dealii_shared_lib_names.insert (p);

  // produce an error if we load deal.II more than once
  if (dealii_shared_lib_names.size() != 1)
  {
    std::ostringstream error;
    error << "........................................................\n"
          << "elASPECT currently links against different versions of the\n"
          << "deal.II library, namely the ones at these locations:\n";
    for (const auto &p : dealii_shared_lib_names)
      error << "  " << p << '\n';
    error << "This can not work.\n\n";

    if (before_loading_shared_libs)
      error << "Since this is happening already before opening additional\n"
            << "shared libraries, this means that something must have gone\n"
            << "wrong when you configured deal.II and/or elASPECT. Please\n"
            << "contact the mailing lists for help.\n";
    else
      error << "Since this is happening after opening additional shared\n"
            << "library plugins, this likely means that you have compiled\n"
            << "elASPECT in release mode and the plugin in debug mode, or the\n"
            << "other way around. Please re-compile the plugin in the same\n"
            << "mode as elASPECT.\n";

    error << "........................................................\n";

    // if not success, then throw an exception: ExcMessage on processor 0,
    // QuietException on the others
    if (dealii::Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      AssertThrow (false, dealii::ExcMessage (error.str()));
    }
    else
      throw elaspect::QuietException();
  }
#else
  // simply mark the argument as read, to avoid compiler warnings
  (void)before_loading_shared_libs;
#endif
}

#endif


// retrieve a list of shared libraries from the parameter file and
// dlopen them so that we can load plugins declared in them
void possibly_load_shared_libs (const std::string &parameters)
{
  using namespace dealii;

  const std::string shared_libs
    = get_last_value_of_parameter(parameters,
                                  "Additional shared libraries");

  if (shared_libs.size() > 0)
  {
#if ELASPECT_USE_SHARED_LIBS==1
    // check up front whether the list of shared libraries is internally
    // consistent or whether we link, for whatever reason, with both the
    // debug and release versions of deal.II
    validate_shared_lib_list (true);

    const std::vector<std::string>
    shared_libs_list = Utilities::split_string_list (shared_libs);

    for (const auto &shared_lib : shared_libs_list)
    {
      if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
        std::cout << "Loading shared library <"
                  << shared_lib
                  << ">" << std::endl;

      void *handle = dlopen (shared_lib.c_str(), RTLD_LAZY);
      AssertThrow (handle != nullptr,
                   ExcMessage (std::string("Could not successfully load shared library <")
                               + shared_lib + ">. The operating system reports "
                               + "that the error is this: <"
                               + dlerror() + ">."));

      // check again whether the list of shared libraries is
      // internally consistent or whether we link with both the
      // debug and release versions of deal.II. this may happen if
      // the plugin was compiled against the debug version of
      // deal.II but elaspect itself against the release version, or
      // the other way around
      validate_shared_lib_list (false);

      // on systems where we can detect that both libdeal_II.so and
      // libdeal_II.g.so is loaded, the test above function above will
      // throw an exception and we will terminate. on the other hand, on
      // systems where we can't detect this we should at least mitigate
      // some of the ill effects -- in particular, make sure that
      // deallog is set to use the desired output depth since otherwise
      // we get lots of output from the linear solvers
      deallog.depth_console(0);
    }

    if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
      std::cout << std::endl;
#else
    if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "You can not load plugins through additional shared libraries " << std::endl
                << "on systems where you link elASPECT as a static executable."
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
    }
    std::exit (1);
#endif
  }
}


std::string
read_parameter_file (const std::string &parameter_file_name)
{
  using namespace dealii;

  std::ifstream parameter_file(parameter_file_name.c_str());
  if (!parameter_file)
  {
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      std::cerr << "Error: Input parameter file <" << parameter_file_name << "> not found."
                << std::endl;

    throw elaspect::QuietException();
    return "";
  }

  std::string input_as_string;
  while (parameter_file)
  {
    std::string line;
    std::getline(parameter_file, line);

    input_as_string += line + '\n';
  }

  return input_as_string;
}


void
parse_parameters (const std::string &input_as_string,
                  dealii::ParameterHandler &prm)
{
  // try reading on processor 0
  bool success = true;
  if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  {
    try
    {
      prm.parse_input_from_string(input_as_string.c_str());
    }
    catch (const dealii::ExceptionBase &e)
    {
      success = false;
      e.print_info(std::cerr);
      std::cerr << std::endl;
    }
  }

  // broadcast the result. we'd like to do this with a bool
  // data type but MPI_C_BOOL is not part of old MPI standards.
  // so, do the broadcast in integers
  {
    int isuccess = (success ? 1 : 0);
    MPI_Bcast (&isuccess, 1, MPI_INT, 0, MPI_COMM_WORLD);
    success = (isuccess == 1);
  }

  // if not success, then throw an exception: ExcMessage on processor 0,
  // QuietException on the others
  if (success == false)
    {
      if (dealii::Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
        {
          AssertThrow(false, dealii::ExcMessage ("Invalid input parameter file."));
        }
      else
        throw elaspect::QuietException();
    }

  // otherwise, processor 0 was ok reading the data, so we can expect the
  // other processors will be ok as well
  if (dealii::Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) != 0)
    {
      prm.parse_input_from_string(input_as_string.c_str());
    }
}


// hook into SIGABRT/SIGFPE and kill off the program
void signal_handler(int signal)
{
  if (signal == SIGABRT)
    {
      std::cerr << "SIGABRT received\n";
    }
  else if (signal == SIGFPE)
    {
      std::cerr << "SIGFPE received\n";
    }
  else
    {
      std::cerr << "Unexpected signal " << signal << " received\n";
    }

  // Kill the program without performing any other cleanup, which would likely
  // lead to a deadlock.
  std::_Exit(EXIT_FAILURE);
}


int main (int argc, char *argv[])
{
  using namespace dealii;

#ifdef DEBUG
#ifdef ELASPECT_USE_FP_EXCEPTIONS
  // Some implementations seem to not initialize the floating point exception
  // bits to zero. Make sure we start from a clean state.
  feclearexcept(FE_DIVBYZERO|FE_INVALID);

  // enable floating point exceptions
  feenableexcept(FE_DIVBYZERO|FE_INVALID);
#endif
#endif

  AssertThrow (argc == 2, 
               ExcMessage("elASPECT accepts only one argument --- the name of the parameter "
                          "file. However, " + Utilities::int_to_string(argc-1) + 
                          " arguments are detected."));

  try
  {
    // Note: we initialize this class inside the try/catch block and not
    // before, so that the destructor of this instance can react if we are
    // currently upwinding the stack if an unhandled exception is being
    // thrown to avoid MPI deadlocks.
    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

    deallog.depth_console(0);

    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      // Output header
      print_elaspect_header(std::cout);
    }
    else
    {
      // We hook into the abort handler on ranks != 0 to avoid an MPI
      // deadlock. The deal.II library will call std::abort() when an
      // Assert is triggered, which can lead to a deadlock because it
      // runs the things that are associated with atexit() which may
      // itself trigger MPI communication. The same happens for other
      // signals we may trigger, such as floating point exceptions
      // (SIGFPE).
      //
      // We work around this by immediately calling _Exit in the
      // signal handler and thus aborting the program without running
      // cleanup functions set via atexit(). This is only necessary on
      // rank != 0 for some reason.
      std::signal(SIGABRT, signal_handler);
      std::signal(SIGFPE, signal_handler);
    }

    const std::string prm_name = argv[1];
    const std::string input_as_string = read_parameter_file(prm_name);
  
    // Determine the dimension we want to work in.
    const unsigned int dim = get_dimension(input_as_string);

    // Do the same with lines potentially indicating shared libs to
    // be loaded. These shared libs could contain additional module
    // instantiations for geometries, etc, that would then be
    // available as part of the possible parameters of the input
    // file, so they need to be loaded before we even start processing
    // the parameter file.
    possibly_load_shared_libs (input_as_string);

    // Now switch between the templates that start the model for 2d or 3d.
    ParameterHandler prm;
    switch (dim)
    {
      case 2:
      {
        elaspect::Simulator<2>::declare_parameters(prm);
        parse_parameters (input_as_string, prm);
      
        elaspect::Simulator<2> simulator (MPI_COMM_WORLD, prm);
        simulator.run();

        break;
      }
      case 3:
      {
        elaspect::Simulator<3>::declare_parameters(prm);
        parse_parameters (input_as_string, prm);
      
        elaspect::Simulator<3> simulator (MPI_COMM_WORLD, prm);
        simulator.run();

        break;
      }
      default:
        AssertThrow((dim >= 2) && (dim <= 3),
                    ExcMessage("elASPECT can only be run in 2d and 3d but a "
                               "different space dimension is given in the parameter file."));
    }
  }
  catch (ExceptionBase &exc)
  {
    // report name of the deal.II exception:
    std::cerr << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception '" << exc.get_exc_name() << "'"
              << " on rank " << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
              << " on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }
  catch (std::exception &exc)
  {
    std::cerr << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception"
              << " on rank " << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
              << " on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }
  catch (elaspect::QuietException &)
  {
    // Quietly treat an exception used on processors other than
    // root when we already know that processor 0 will generate
    // an exception. We do this to avoid creating too much
    // (duplicate) screen output.

    // Sleep a few seconds before aborting. This allows text output from
    // other ranks to be printed before the MPI implementation might kill
    // the computation.
    std::this_thread::sleep_for(std::chrono::seconds(5));

    MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }
  catch (...)
  {
    std::cerr << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }

  return 0;
}
