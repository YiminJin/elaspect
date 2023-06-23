#include <elaspect/parameters.h>
#include <elaspect/simulator.h>

#include <boost/lexical_cast.hpp>

namespace elaspect
{
  template <int dim>
  Parameters<dim>::Parameters (ParameterHandler &prm,
                               MPI_Comm mpi_communicator)
  {
    parse_parameters (prm, mpi_communicator);
  }


  template <int dim>
  void
  Parameters<dim>::
  declare_parameters (ParameterHandler &prm)
  {
    prm.declare_entry ("Dimension", "2",
                       Patterns::Integer (2,3),
                       "The number of space dimensions you want to run this program in. "
                       "elASPECT can run in 2 and 3 space dimensions.");

    prm.declare_entry ("Additional shared libraries", "",
                       Patterns::List (Patterns::FileName()),
                       "A list of names of additional shared libraries that should be loaded "
                       "upon starting up the program. The names of these files can contain absolute "
                       "or relative paths (relative to the directory in which you call elASPECT). "
                       "In fact, file names that do not contain any directory "
                       "information (i.e., only the name of a file such as <myplugin.so> "
                       "will not be found if they are not located in one of the directories "
                       "listed in the \\texttt{LD_LIBRARY_PATH} environment variable. In order "
                       "to load a library in the current directory, use <./myplugin.so> "
                       "instead."
                       "\n\n"
                       "The typical use of this parameter is so that you can implement "
                       "additional plugins in your own directories, rather than in the elASPECT "
                       "source directories. You can then simply compile these plugins into a "
                       "shared library without having to re-compile all of elASPECT. See the "
                       "section of the manual discussing writing extensions for more "
                       "information on how to compile additional files into a shared "
                       "library.");

    prm.declare_entry ("Output directory", "output",
                       Patterns::DirectoryName(),
                       "The name of the directory into which all output files should be "
                       "placed. This may be an absolute or a relative path.");

    prm.declare_entry ("Use years in output instead of seconds", "true",
                       Patterns::Bool (),
                       "");

    prm.declare_entry ("Start time", "0.",
                       Patterns::Double (),
                       "The start time of the simulation. Units: Years if the "
                       "'Use years in output instead of seconds' parameter is set; "
                       "seconds otherwise.");

    prm.declare_entry ("Timing output frequency", "100",
                       Patterns::Integer(0),
                       "How frequently in timesteps to output timing information. This is "
                       "generally adjusted only for debugging and timing purposes. If the "
                       "value is set to zero it will also output timing information at the "
                       "initiation timesteps.");

    prm.declare_entry ("Use ALE method", "false",
                       Patterns::Bool(),
                       "");

    prm.declare_entry ("Include heat transport", "false",
                       Patterns::Bool(),
                       "");

    prm.declare_entry ("CFL number", "1.0",
                       Patterns::Double(0),
                       "");

    prm.enter_subsection("Compositional fields");
    {
      prm.declare_entry ("Number of fields", "0",
                         Patterns::Integer(0),
                         "The number of fields that will be advected with the mesh, excluding displacement, "
                         "temperature and stress.");
      prm.declare_entry ("Names of fields", "",
                         Patterns::List(Patterns::Anything()),
                         "A user-defined name for each of the compositional fields requested.");
    }
    prm.leave_subsection();

    prm.enter_subsection("Solver parameters");
    {
      prm.enter_subsection("Mechanical system");
      {
        prm.declare_entry ("Maximum nonlinear iterations", "10",
                           Patterns::Integer(1),
                           "The maximal number of nonlinear iterations to be performed. "
                           "This parameter is active only when plasticity is included "
                           "in the constitutive relations.");

        prm.declare_entry ("Nonlinear solver tolerance", "1e-8",
                           Patterns::Double(0, 1),
                           "A relative tolerance up to which the nonlinear solver will iterate. "
                           "This parameter is active only when plasticity is included in the "
                           "constitutive relations.");

        prm.declare_entry ("Maximum linear iterations", "5000",
                           Patterns::Integer(1),
                           "");

        prm.declare_entry ("Linear solver tolerance", "1e-8",
                           Patterns::Double(0,1),
                           "");

        prm.declare_entry ("Linear solver scheme", "CG",
                           Patterns::Selection("CG|Bicgstab|GMRES|MUMPS"),
                           "");

        prm.declare_entry ("GMRES solver restart length", "200",
                           Patterns::Integer(1),
                           "");

        prm.declare_entry ("AMG aggregation threshold", "0.01",
                           Patterns::Double(0, 1),
                           "");

        prm.declare_entry ("Maximum line search steps", "3",
                           Patterns::Integer(0),
                           "");

        prm.declare_entry ("Line search beta", "1e-4",
                           Patterns::Double(0, 0.5),
                           "");

        prm.declare_entry ("Line search delta", "0.5",
                           Patterns::Double(0, 1),
                           "");

        prm.declare_entry ("Enforce convergence", "false",
                           Patterns::Bool(),
                           "");
      }
      prm.leave_subsection();

      prm.enter_subsection("Thermal system");
      {
        prm.declare_entry ("Linear solver tolerance", "1e-12",
                           Patterns::Double(0,1),
                           "");

        prm.declare_entry ("Maximum linear iterations", "1000",
                           Patterns::Integer(0),
                           "");
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();

    prm.enter_subsection ("Nullspace removal");
    {
      prm.declare_entry ("Remove nullspace", "",
                         Patterns::MultipleSelection("net rotation|angular momentum|"
                                                     "net translation|linear momentum|"
                                                     "net x translation|net y translation|net z translation|"
                                                     "linear x momentum|linear y momentum|linear z momentum"),
                         "Choose none, one or several from "
                         "\n\n"
                         "\\begin{itemize} \\item net rotation \\item angular momentum \\item net translation "
                         "\\item linear momentum \\item net x translation \\item net y translation "
                         "\\item net z translation \\item linear x momentum \\item linear y momentum "
                         "\\item linear z momentum \\end{itemize}"
                         "\n\n"
                         "These are a selection of operations to remove certain parts of the nullspace from "
                         "the displacement after solving. For some geometries and certain boundary conditions "
                         "the displacement field is not uniquely determined but contains free translations "
                         "and/or rotations. Depending on what you specify here, these non-determined "
                         "modes will be removed from the displacement field at the end of the Stokes solve step.\n"
                         "\n\n"
                         "The ``angular momentum'' option removes a rotation such that the net angular momentum "
                         "is zero. The ``linear * momentum'' options remove translations such that the net "
                         "momentum in the relevant direction is zero.  The ``net rotation'' option removes the "
                         "net rotation of the whole domain, and the ``net * translation'' options remove the "
                         "net translations in the relevant directions.  For most problems there should not be a "
                         "significant difference between the momentum and rotation/translation versions of "
                         "nullspace removal, although the momentum versions are more physically motivated. "
                         "They are equivalent for constant density simulations, and approximately equivalent "
                         "when the density variations are small."
                         "\n\n"
                         "Note that while more than one operation can be selected it only makes sense to "
                         "pick one rotational and one translational operation.");
    }
    prm.leave_subsection();

    prm.enter_subsection("Mesh refinement");
    {
      prm.declare_entry ("Initial global refinement", "2",
                         Patterns::Integer (0),
                         "The number of global refinement steps performed on "
                         "the initial coarse mesh, before the problem is first "
                         "solved there.");

      prm.declare_entry ("Initial adaptive refinement", "0",
                         Patterns::Integer (0),
                         "The number of adaptive refinement steps performed after "
                         "initial global refinement but while still within the first "
                         "time step. These refinement steps (n) are added to the value "
                         "for initial global refinement (m) so that the final mesh has "
                         "cells that are at most on refinement level $n+m$.");

      prm.declare_entry ("Time steps between mesh refinement", "10",
                         Patterns::Integer (0),
                         "The number of time steps after which the mesh is to be "
                         "adapted again based on computed error indicators. If 0 "
                         "then the mesh will never be changed.");

      prm.declare_entry ("Refinement fraction", "0.3",
                         Patterns::Double(0., 1.),
                         "Cells are sorted from largest to smallest by their total error "
                         "(determined by the Strategy). Then the cells with the largest "
                         "error (top of this sorted list) that account for given fraction "
                         "of the error are refined.");

      prm.declare_entry ("Coarsening fraction", "0.05",
                         Patterns::Double(0., 1.),
                         "Cells are sorted from largest to smallest by their total error "
                         "(determined by the Strategy). Then the cells with the smallest "
                         "error (bottom of this sorted list) that account for the given fraction "
                         "of the error are coarsened.");

      prm.declare_entry ("Minimum refinement level", "0",
                         Patterns::Integer (0),
                         "The minimum refinement level each cell should have, "
                         "and that can not be exceeded by coarsening. "
                         "Should not be higher than the 'Initial global refinement' "
                         "parameter.");

      prm.declare_entry ("Adapt by fraction of cells", "false",
                         Patterns::Bool(),
                         "Use fraction of the total number of cells instead of "
                         "fraction of the total error as the limit for refinement "
                         "and coarsening.");

      prm.declare_entry ("Skip solvers on initial refinement", "false",
                         Patterns::Bool (),
                         "Whether or not solvers should be executed during the initial "
                         "adaptive refinement cycles that are run at the start of the "
                         "simulation.");

      prm.declare_entry ("Run postprocessors on initial refinement", "false",
                         Patterns::Bool (),
                         "Whether or not the postprocessors should be executed after "
                         "each of the initial adaptive refinement cycles that are run at "
                         "the start of the simulation. This is useful for "
                         "plotting/analyzing how the mesh refinement parameters are "
                         "working for a particular model.");
    }
    prm.leave_subsection();

    prm.enter_subsection("Discretization");
    {
      prm.declare_entry("Displacement polynomial degree", "2",
                        Patterns::Integer(1),
                        "");

      prm.declare_entry("Fluid pressure polynomial degree", "1",
                        Patterns::Integer(1),
                        "");

      prm.declare_entry("Temperature polynomial degree", "2",
                        Patterns::Integer(1),
                        "");

      prm.declare_entry("Composition polynomial degree", "2",
                        Patterns::Integer(1),
                        "");

      prm.declare_entry("Number of Gaussian points per dimension", "3",
                        Patterns::Integer(1),
                        "");
    }
    prm.leave_subsection();

    prm.enter_subsection("Stabilization parameters");
    {
      prm.declare_entry("Apply BP limiter to compositional fields", "false",
                        Patterns::Bool(),
                        "Whether to apply the bound preserving limiter to compositional fields "
                        "as a correction after advection when ALE method is applied. This limiter "
                        "keeps the discontinuous solution in the range given by "
                        "'Global composition maximum' and 'Global composition minimum'.");

      prm.declare_entry("Global composition maximum",
                        boost::lexical_cast<std::string>(std::numeric_limits<double>::max()),
                        Patterns::List(Patterns::Double()),
                        "The maximum global composition values that will be used in the bound "
                        "preserving limiter for composition advection. The number of input values "
                        "separated by ',' has to be one or the same as the number of the "
                        "compositional fields. When only one value is supplied, the same value is "
                        "assumed for all compositional fields.");

      prm.declare_entry("Global composition minimum",
                        boost::lexical_cast<std::string>(std::numeric_limits<double>::lowest()),
                        Patterns::List(Patterns::Double()),
                        "The minimum global composition values that will be used in the bound "
                        "preserving limiter for composition advection. The number of input values "
                        "separated by ',' has to be one or the same as the number of the "
                        "compositional fields. When only one value is supplied, the same value is "
                        "assumed for all compositional fields.");

      prm.declare_entry("Apply WENO limiter to physical fields", "false",
                        Patterns::Bool(),
                        "Whether to apply WENO limiter to physical fields in quadrature point data "
                        "as a correction after advection when ALE method is applied. This limiter "
                        "replace the DG polynomials in troubled cells with reconstructed "
                        "polynomials that keep the original cell averages and are less oscillatory. "
                        "This limiter can only be applied when the polynomial degree of physical "
                        "fields in quadrature point data is 1.");

      prm.declare_entry("KXRCF indicator threshold", "1.0",
                        Patterns::Double(0),
                        "The threshold $C_k$ in KXRCF indicator. In the WENO limiting procedure, "
                        "a cell is marked as troubled cell if its indicating quantity for KXRCF "
                        "is greater than $C_k$.");
    }
    prm.leave_subsection();

    prm.enter_subsection("Postprocess");
    {
      prm.declare_entry ("Run postprocessors on nonlinear iterations", "false",
                         Patterns::Bool (),
                         "Whether or not the postprocessors should be executed after "
                         "each of the nonlinear iterations done within one time step. "
                         "As this is mainly an option for the purposes of debugging, "
                         "it is not supported when the 'Time between graphical output' "
                         "is larger than zero, or when the postprocessor is not intended "
                         "to be run more than once per timestep.");
    }
    prm.leave_subsection();

    prm.enter_subsection ("Boundary traction model");
    {
      prm.declare_entry ("Prescribed traction boundary indicators", "",
                         Patterns::Map (Patterns::Anything(),
                                        Patterns::Selection(BoundaryTraction::get_names<dim>())),
                         "A comma separated list denoting those boundaries "
                         "on which a traction force is prescribed, i.e., where "
                         "known external forces act, resulting in an unknown velocity. This is "
                         "often used to model ``open'' boundaries where we only know the pressure. "
                         "This pressure then produces a force that is normal to the boundary and "
                         "proportional to the pressure."
                         "\n\n"
                         "The format of valid entries for this parameter is that of a map "
                         "given as ``key1 [selector]: value1, key2 [selector]: value2, key3: value3, ...'' where "
                         "each key must be a valid boundary indicator (which is either an "
                         "integer or the symbolic name the geometry model in use may have "
                         "provided for this part of the boundary) "
                         "and each value must be one of the currently implemented boundary "
                         "traction models. ``selector'' is an optional string given as a subset "
                         "of the letters `xyz' that allows you to apply the boundary conditions "
                         "only to the components listed. As an example, '1 y: function' applies "
                         "the type `function' to the y component on boundary 1. Without a selector "
                         "it will affect all components of the traction.");
    }
    prm.leave_subsection();

    prm.enter_subsection ("Boundary heat flux model");
    {
      prm.declare_entry ("Fixed heat flux boundary indicators", "",
                         Patterns::List (Patterns::Anything()),
                         "A comma separated list of names denoting those boundaries "
                         "on which the heat flux is fixed and described by the "
                         "boundary heat flux object selected in the 'Model name' parameter. "
                         "All boundary indicators used by the geometry but not explicitly "
                         "listed here or in the list of 'Fixed temperature boundary indicators' "
                         "in the 'Boundary temperature model' will end up with no-flux "
                         "(insulating) boundary conditions."
                         "\n\n"
                         "The names of the boundaries listed here can either be "
                         "numbers (in which case they correspond to the numerical "
                         "boundary indicators assigned by the geometry object), or they "
                         "can correspond to any of the symbolic names the geometry object "
                         "may have provided for each part of the boundary. You may want "
                         "to compare this with the documentation of the geometry model you "
                         "use in your model."
                         "\n\n"
                         "This parameter only describes which boundaries have a fixed "
                         "heat flux, but not what heat flux should hold on these "
                         "boundaries. The latter piece of information needs to be "
                         "implemented in a plugin in the BoundaryHeatFlux "
                         "group, unless an existing implementation in this group "
                         "already provides what you want.");
    }
    prm.leave_subsection();


    prm.enter_subsection ("Material model");
    {
      prm.declare_entry ("Constitutive relation", "elasticity",
                         Patterns::List(Patterns::Selection("elasticity|viscosity|plasticity|thermal expansion|pore fluid")),
                         "");

      prm.declare_entry ("Material averaging", "none",
                         Patterns::Selection(MaterialModel::MaterialAveraging::
                                             get_averaging_operation_names()),
                         "Whether or not (and in the first case, how) to do any averaging of "
                         "material model output data when constructing the linear systems "
                         "for velocity/pressure, temperature, and compositions in each "
                         "time step, as well as their corresponding preconditioners."
                         "\n\n"
                         "Possible choices: " + MaterialModel::MaterialAveraging::
                         get_averaging_operation_names()
                         +
                         "\n\n"
                         "The process of averaging, and where it may be used, is "
                         "discussed in more detail in "
                         "Section~\\ref{sec:sinker-with-averaging}."
                         "\n\n"
                         "More averaging schemes are available in the averaging material "
                         "model. This material model is a ``compositing material model'' "
                         "which can be used in combination with other material models.");
    }
    prm.leave_subsection ();
  }


  template <int dim>
  void
  Parameters<dim>::
  parse_parameters (ParameterHandler &prm,
                    const MPI_Comm mpi_communicator)
  {
    // first, make sure that the ParameterHandler parser agrees
    // with the code in main() about the meaning of the "Dimension"
    // parameter
    AssertThrow (prm.get_integer("Dimension") == dim,
                 ExcInternalError());

    output_directory = prm.get("Output directory");
    if (output_directory.size() == 0)
      output_directory = "./";
    else if (output_directory[output_directory.size()-1] != '/')
      output_directory += "/";

    Utilities::create_directory (output_directory,
                                 mpi_communicator,
                                 false);

    convert_to_years        = prm.get_bool ("Use years in output instead of seconds");
    start_time              = prm.get_double ("Start time");
    timing_output_frequency = prm.get_integer ("Timing output frequency");
    use_ALE_method          = prm.get_bool ("Use ALE method");
    include_heat_transport  = prm.get_bool ("Include heat transport");
    CFL_number              = prm.get_double ("CFL number");

    prm.enter_subsection("Compositional fields");
    {
      n_compositional_fields  = prm.get_integer ("Number of fields");

      names_of_compositional_fields = Utilities::split_string_list(prm.get("Names of fields"));
      AssertThrow ((names_of_compositional_fields.size() == 0) ||
                   (names_of_compositional_fields.size() == n_compositional_fields),
                   ExcMessage("The length of the list of names for the compositional "
                              "fields needs to either be empty or have length equal to "
                              "the number of compositional fields."));

      // check that the names use only allowed characters, are not empty strings and are unique
      for (unsigned int i = 0; i < names_of_compositional_fields.size(); ++i)
      {
        AssertThrow (names_of_compositional_fields[i].find_first_not_of("abcdefghijklmnopqrstuvwxyz"
                                                                        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                                                        "0123456789_") == std::string::npos,
                     ExcMessage("Invalid character in field " + names_of_compositional_fields[i] + ". "
                                "Names of compositional fields should consist of a "
                                "combination of letters, numbers and underscores."));
        AssertThrow (names_of_compositional_fields[i].size() > 0,
                     ExcMessage("Invalid name of field " + names_of_compositional_fields[i] + ". "
                                "Names of compositional fields need to be non-empty."));

        for (unsigned int j = 0; j < i; ++j)
          AssertThrow (names_of_compositional_fields[i] != names_of_compositional_fields[j],
                       ExcMessage("Names of compositional fields have to be unique! " + names_of_compositional_fields[i] +
                                  " is used more than once."));
      }

      // default names if list is empty
      if (names_of_compositional_fields.size() == 0)
        for (unsigned int i = 0; i < n_compositional_fields; ++i)
          names_of_compositional_fields.push_back("C_" + Utilities::int_to_string(i+1));
    }
    prm.leave_subsection();

    prm.enter_subsection("Solver parameters");
    {

      prm.enter_subsection("Mechanical system");
      {
        max_nonlinear_iterations = prm.get_integer ("Maximum nonlinear iterations");
        nonlinear_tolerance = prm.get_double ("Nonlinear solver tolerance");

        const std::string x_linear_solver = prm.get("Linear solver scheme");
        if (x_linear_solver == "CG")
          mechanical_system_linear_solver = LinearSolver::CG;
        else if (x_linear_solver == "Bicgstab")
          mechanical_system_linear_solver = LinearSolver::BiCGStab;
        else if (x_linear_solver == "GMRES")
          mechanical_system_linear_solver = LinearSolver::GMRES;
        else if (x_linear_solver == "MUMPS")
          mechanical_system_linear_solver = LinearSolver::MUMPS;
        else
          AssertThrow (false, ExcNotImplemented());

        mechanical_system_solver_tolerance       = prm.get_double("Linear solver tolerance");
        mechanical_system_max_linear_iterations  = prm.get_integer("Maximum linear iterations");
        mechanical_system_gmres_restart_length   = prm.get_integer("GMRES solver restart length");
        mechanical_system_amg_aggregation_threshold = prm.get_double ("AMG aggregation threshold");
        
        max_line_search_steps = prm.get_integer("Maximum line search steps");
        line_search_beta      = prm.get_double("Line search beta");
        line_search_delta     = prm.get_double("Line search delta");

        enforce_convergence_for_mechanical_system = prm.get_bool("Enforce convergence");
      }
      prm.leave_subsection();

      prm.enter_subsection("Thermal system");
      {
        thermo_system_solver_tolerance      = prm.get_double("Linear solver tolerance");
        thermo_system_max_linear_iterations = prm.get_integer("Maximum linear iterations");
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();

    prm.enter_subsection ("Nullspace removal");
    {
      nullspace_removal = NullspaceRemoval::none;
      std::vector<std::string> nullspace_names =
        Utilities::split_string_list(prm.get("Remove nullspace"));
      AssertThrow(Utilities::has_unique_entries(nullspace_names),
                  ExcMessage("The list of strings for the parameter "
                             "'Nullspace removal/Remove nullspace' contains entries more than once. "
                             "This is not allowed. Please check your parameter file."));

      for (unsigned int i=0; i<nullspace_names.size(); ++i)
      {
        if (nullspace_names[i]=="net rotation")
          nullspace_removal = typename NullspaceRemoval::Kind(
                                nullspace_removal | NullspaceRemoval::net_rotation);
        else if (nullspace_names[i]=="angular momentum")
          nullspace_removal = typename NullspaceRemoval::Kind(
                                nullspace_removal | NullspaceRemoval::angular_momentum);
        else if (nullspace_names[i]=="net translation")
          nullspace_removal = typename NullspaceRemoval::Kind(
                                nullspace_removal | NullspaceRemoval::net_translation_x |
                                NullspaceRemoval::net_translation_y | ( dim == 3 ?
                                                                        NullspaceRemoval::net_translation_z : 0) );
        else if (nullspace_names[i]=="net x translation")
          nullspace_removal = typename NullspaceRemoval::Kind(
                                nullspace_removal | NullspaceRemoval::net_translation_x);
        else if (nullspace_names[i]=="net y translation")
          nullspace_removal = typename NullspaceRemoval::Kind(
                                nullspace_removal | NullspaceRemoval::net_translation_y);
        else if (nullspace_names[i]=="net z translation")
          nullspace_removal = typename NullspaceRemoval::Kind(
                                nullspace_removal | NullspaceRemoval::net_translation_z);
        else if (nullspace_names[i]=="linear x momentum")
          nullspace_removal = typename NullspaceRemoval::Kind(
                                nullspace_removal | NullspaceRemoval::linear_momentum_x);
        else if (nullspace_names[i]=="linear y momentum")
          nullspace_removal = typename NullspaceRemoval::Kind(
                                nullspace_removal | NullspaceRemoval::linear_momentum_y);
        else if (nullspace_names[i]=="linear z momentum")
          nullspace_removal = typename NullspaceRemoval::Kind(
                                nullspace_removal | NullspaceRemoval::linear_momentum_z);
        else if (nullspace_names[i]=="linear momentum")
          nullspace_removal = typename NullspaceRemoval::Kind(
                                nullspace_removal | NullspaceRemoval::linear_momentum_x |
                                NullspaceRemoval::linear_momentum_y | ( dim == 3 ?
                                                                        NullspaceRemoval::linear_momentum_z : 0) );
        else
          AssertThrow(false, ExcInternalError());
      }
    }
    prm.leave_subsection ();

    prm.enter_subsection("Mesh refinement");
    {
      initial_global_refinement          = prm.get_integer("Initial global refinement");
      initial_adaptive_refinement        = prm.get_integer("Initial adaptive refinement");
      refinement_fraction                = prm.get_double("Refinement fraction");
      coarsening_fraction                = prm.get_double("Coarsening fraction");
      min_grid_level                     = prm.get_integer ("Minimum refinement level");
      adapt_by_fraction_of_cells         = prm.get_bool ("Adapt by fraction of cells");
      skip_solvers_on_initial_refinement = prm.get_bool("Skip solvers on initial refinement");
      adaptive_refinement_interval       = prm.get_integer ("Time steps between mesh refinement");
      run_postprocessors_on_initial_refinement = prm.get_bool("Run postprocessors on initial refinement");
    }
    prm.leave_subsection();

    prm.enter_subsection("Discretization");
    {
      displacement_degree   = prm.get_integer("Displacement polynomial degree");
      fluid_pressure_degree = prm.get_integer("Fluid pressure polynomial degree");
      temperature_degree    = prm.get_integer("Temperature polynomial degree");
      composition_degree    = prm.get_integer("Composition polynomial degree"); 
      n_gaussian_points     = prm.get_integer("Number of Gaussian points per dimension");
    }
    prm.leave_subsection();

    prm.enter_subsection("Stabilization parameters");
    {
      apply_BP_limiter_to_compositional_fields = prm.get_bool("Apply BP limiter to compositional fields");
      apply_WENO_limiter_to_physical_fields    = prm.get_bool("Apply WENO limiter to physical fields");
      KXRCF_indicator_threshold                = prm.get_double("KXRCF indicator threshold");

      composition_max_preset = Utilities::possibly_extend_from_1_to_N(
        Utilities::string_to_double(Utilities::split_string_list(prm.get("Global composition maximum"))),
        n_compositional_fields,
        "Global composition maximum");

      composition_min_preset = Utilities::possibly_extend_from_1_to_N(
        Utilities::string_to_double(Utilities::split_string_list(prm.get("Global composition minimum"))),
        n_compositional_fields,
        "Global composition minimum");
    }
    prm.leave_subsection();

    prm.enter_subsection ("Postprocess");
    {
      run_postprocessors_on_nonlinear_iterations = prm.get_bool("Run postprocessors on nonlinear iterations");
    }
    prm.leave_subsection ();

    prm.enter_subsection ("Material model");
    {
      std::vector<std::string> constitutive_relation_components = 
        Utilities::split_string_list(prm.get("Constitutive relation"));

      AssertThrow (std::find(constitutive_relation_components.begin(),
                             constitutive_relation_components.end(),
                             std::string("elasticity"))
                   != constitutive_relation_components.end(),
                   ExcMessage("Elasticity must be included in constitutive relations"));

      constitutive_relation = ConstitutiveRelation::elasticity;
      for (const auto &component : constitutive_relation_components)
      {
        if (component == "elasticity")
          continue;
        else if (component == "viscosity")
          constitutive_relation |= ConstitutiveRelation::viscosity;
        else if (component == "plasticity")
          constitutive_relation |= ConstitutiveRelation::plasticity;
        else if (component == "thermal expansion")
          constitutive_relation |= ConstitutiveRelation::thermal_expansion;
        else if (component == "pore fluid")
          constitutive_relation |= ConstitutiveRelation::pore_fluid;
        else
          AssertThrow (false, ExcNotImplemented());
      }

      material_averaging
        = MaterialModel::MaterialAveraging::parse_averaging_operation_name
          (prm.get("Material averaging"));
    }
    prm.leave_subsection ();
  }


  template <int dim>
  void
  Parameters<dim>::
  parse_geometry_dependent_parameters (ParameterHandler &prm,
                                       const GeometryModel::Interface<dim> &geometry_model)
  {
    prm.enter_subsection("Boundary traction model");
    {
      const std::vector<std::string> x_prescribed_traction_boundary_indicators
        = Utilities::split_string_list
          (prm.get ("Prescribed traction boundary indicators"));
      for (const auto &p : x_prescribed_traction_boundary_indicators)
      {
        // each entry has the format (white space is optional):
        // <id> [x][y][z] : <value (might have spaces)>
        //
        // first tease apart the two halves
        const std::vector<std::string> split_parts = Utilities::split_string_list (p, ':');
        AssertThrow (split_parts.size() == 2,
                     ExcMessage ("The format for prescribed traction boundary indicators "
                                 "requires that each entry has the form `"
                                 "<id> [x][y][z] : <value>', but there does not "
                                 "appear to be a colon in the entry <"
                                 + p
                                 + ">."));

        // the easy part: get the value
        const std::string value = split_parts[1];

        // now for the rest. since we don't know whether there is a
        // component selector, start reading at the end and subtracting
        // letters x, y and z
        std::string key_and_comp = split_parts[0];
        std::string comp;
        while ((key_and_comp.size()>0) &&
               ((key_and_comp[key_and_comp.size()-1] == 'x')
                ||
                (key_and_comp[key_and_comp.size()-1] == 'y')
                ||
                ((key_and_comp[key_and_comp.size()-1] == 'z') && (dim==3))))
        {
          comp += key_and_comp[key_and_comp.size()-1];
          key_and_comp.erase (--key_and_comp.end());
        }

        // we've stopped reading component selectors now. there are three
        // possibilities:
        // - no characters are left. this means that key_and_comp only
        //   consisted of a single word that only consisted of 'x', 'y'
        //   and 'z's. then this would have been a mistake to classify
        //   as a component selector, and we better undo it
        // - the last character of key_and_comp is not a whitespace. this
        //   means that the last word in key_and_comp ended in an 'x', 'y'
        //   or 'z', but this was not meant to be a component selector.
        //   in that case, put these characters back.
        // - otherwise, we split successfully. eat spaces that may be at
        //   the end of key_and_comp to get key
        if (key_and_comp.size() == 0)
          key_and_comp.swap (comp);
        else if (key_and_comp[key_and_comp.size()-1] != ' ')
        {
          key_and_comp += comp;
          comp = "";
        }
        else
        {
          while ((key_and_comp.size()>0) && (key_and_comp[key_and_comp.size()-1] == ' '))
            key_and_comp.erase (--key_and_comp.end());
        }

        // finally, try to translate the key into a boundary_id. then
        // make sure we haven't seen it yet
        types::boundary_id boundary_id;
        try
        {
          boundary_id = geometry_model.translate_symbolic_boundary_name_to_id(key_and_comp);
        }
        catch (const std::string &error)
        {
          AssertThrow (false, ExcMessage ("While parsing the entry <Boundary traction model/Prescribed "
                                          "traction indicators>, there was an error. Specifically, "
                                          "the conversion function complained as follows:\n\n"
                                          + error));
        }

        AssertThrow (prescribed_traction_boundary_indicators.find(boundary_id)
                     == prescribed_traction_boundary_indicators.end(),
                     ExcMessage ("Boundary indicator <" + Utilities::int_to_string(boundary_id) +
                                 "> appears more than once in the list of indicators "
                                 "for nonzero traction boundaries."));

        // finally, put it into the list
        prescribed_traction_boundary_indicators[boundary_id] =
          std::pair<std::string,std::string>(comp,value);
      }
    }
    prm.leave_subsection ();

    prm.enter_subsection ("Boundary heat flux model");
    {
      try
      {
        const std::vector<types::boundary_id> x_fixed_heat_flux_boundary_indicators
          = geometry_model.translate_symbolic_boundary_names_to_ids(Utilities::split_string_list
                                                                    (prm.get ("Fixed heat flux boundary indicators")));
        fixed_heat_flux_boundary_indicators
          = std::set<types::boundary_id> (x_fixed_heat_flux_boundary_indicators.begin(),
                                          x_fixed_heat_flux_boundary_indicators.end());
      }
      catch (const std::string &error)
      {
        AssertThrow (false, ExcMessage ("While parsing the entry <Boundary heat flux model/Fixed heat flux "
                                        "boundary indicators>, there was an error. Specifically, "
                                        "the conversion function complained as follows:\n\n"
                                        + error));
      }
    }
    prm.leave_subsection ();
  }


  template <int dim>
  void
  Simulator<dim>::declare_parameters(ParameterHandler &prm)
  {
    Parameters<dim>::declare_parameters(prm);
    Postprocess::Manager<dim>::declare_parameters(prm);
    MeshDeformationHandler<dim>::declare_parameters(prm);
    MeshRefinement::Manager<dim>::declare_parameters(prm);
    GeometryModel::declare_parameters<dim>(prm);
    InitialTopography::declare_parameters<dim>(prm);
    HeatingModel::Manager<dim>::declare_parameters(prm);
    MaterialHandler<dim>::declare_parameters(prm);
    GravityModel::declare_parameters<dim>(prm);
    TimeStepping::Manager<dim>::declare_parameters(prm);
    InitialComposition::Manager<dim>::declare_parameters(prm);
    InitialTemperature::Manager<dim>::declare_parameters(prm);
    BoundaryComposition::Manager<dim>::declare_parameters(prm);
    BoundaryTemperature::Manager<dim>::declare_parameters(prm);
    BoundaryVelocity::Manager<dim>::declare_parameters(prm);
    BoundaryTraction::declare_parameters<dim>(prm);
    BoundaryHeatFlux::declare_parameters<dim>(prm);
  }
}


// explicit instantiations
namespace elaspect
{
#define INSTANTIATE(dim) \
  template struct Parameters<dim>; \
  template void Simulator<dim>::declare_parameters (ParameterHandler &prm);

  ELASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
