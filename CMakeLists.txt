CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

MESSAGE(STATUS "====================================================")
MESSAGE(STATUS "============== Configuring elASPECT=================")
MESSAGE(STATUS "====================================================")

SET(TARGET "elaspect")

FILE(GLOB_RECURSE TARGET_SRC "source/*.cc" "include/*.h")

# Set up include directories.
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include include)

SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

LIST(APPEND CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/cmake/modules
  )

FIND_PACKAGE(deal.II 9.2.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../  $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n*** Could not find a suitably recent version of deal.II. ***\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake "
    "or set an environment variable \"DEAL_II_DIR\" that contains a path to a "
    "sufficiently recent version of deal.II."
    )
ENDIF()

MESSAGE(STATUS "Found deal.II version ${DEAL_II_PACKAGE_VERSION} at '${deal.II_DIR}'")

IF(NOT DEAL_II_WITH_CXX11)
  MESSAGE(FATAL_ERROR "\n*** elASPECT requires C++11 but your version of deal.II is not "
    "configured with it. This likely means that your compiler is too old. Check "
    "DEAL_II_WITH_CXX11 in deal.II.")
ENDIF()

SET(_DEALII_GOOD ON)

IF(NOT DEAL_II_WITH_P4EST OR NOT DEAL_II_WITH_TRILINOS)
  MESSAGE(FATAL_ERROR
    "\nelASPECT requires a deal.II installation built with support for Trilinos and p4est but one or both of these appears to be missing!\n"
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET} CXX)

# Pass down the source directory to the sources.
FOREACH(_source_file ${TARGET_SRC})
  SET_PROPERTY(SOURCE ${_source_file}
    APPEND PROPERTY COMPILE_DEFINITIONS ELASPECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
ENDFOREACH()

# load in version info and export it
FILE(STRINGS "${CMAKE_SOURCE_DIR}/VERSION" ELASPECT_PACKAGE_VERSION LIMIT_COUNT 1)

INCLUDE(CMakePackageConfigHelpers)
WRITE_BASIC_PACKAGE_VERSION_FILE(
  "${CMAKE_BINARY_DIR}/ElaspectConfigVersion.cmake"
  VERSION ${ELASPECT_PACKAGE_VERSION}
  COMPATIBILITY AnyNewerVersion
)

# Configure a cmake fragment that plugins can use to
# set up compiler flags, include paths, etc to compile an
# elASPECT plugin.
SET(CONFIG_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/include")
SET(CONFIG_DIR "${CMAKE_BINARY_DIR}")
CONFIGURE_FILE(
  ${CMAKE_SOURCE_DIR}/cmake/ElaspectConfig.cmake.in
  ${CMAKE_BINARY_DIR}/ElaspectConfig.cmake
  @ONLY
)

# Config for the install dir:
SET(CONFIG_INCLUDE_DIRS "${CMAKE_INSTALL_PREFIX}/include")
SET(CONFIG_DIR "${CMAKE_INSTALL_PREFIX}/bin")
CONFIGURE_FILE(
  ${CMAKE_SOURCE_DIR}/cmake/ElaspectConfig.cmake.in
  ${CMAKE_BINARY_DIR}/forinstall/ElaspectConfig.cmake
  @ONLY
  )

# Provide "release" and "debug" targets to switch compile mode
IF(${DEAL_II_BUILD_TYPE} MATCHES "DebugRelease")
  ADD_CUSTOM_TARGET(release
    COMMAND ${CMAKE_COMMAND} -D CMAKE_BUILD_TYPE=Release .
    COMMAND ${CMAKE_COMMAND} -E echo "***"
    COMMAND ${CMAKE_COMMAND} -E echo "*** Switched to Release mode. Now recompile with: ${_make_command}"
    COMMAND ${CMAKE_COMMAND} -E echo "***"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    VERBATIM
    COMMENT "switching to RELEASE mode..."
    )
  ADD_CUSTOM_TARGET(debug
    COMMAND ${CMAKE_COMMAND} -D CMAKE_BUILD_TYPE=Debug .
    COMMAND ${CMAKE_COMMAND} -E echo "***"
    COMMAND ${CMAKE_COMMAND} -E echo "*** Switched to Debug mode. Now recompile with: ${_make_command}"
    COMMAND ${CMAKE_COMMAND} -E echo "***"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    VERBATIM
    COMMENT "switching to DEBUG mode..."
    )
ENDIF()

# Provide a "distclean" target (like it is done in deal.II):
ADD_CUSTOM_TARGET(distclean
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
  COMMAND ${CMAKE_COMMAND} -E remove_directory CMakeFiles
  COMMAND ${CMAKE_COMMAND} -E remove
    CMakeCache.txt cmake_install.cmake Makefile
    build.ninja rules.ninja .ninja_deps .ninja_log
  COMMENT "distclean invoked"
  )

FILE(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake
"MESSAGE(
\"###
#
#  Project ${TARGET} set up with  ${DEAL_II_PACKAGE_NAME}-${DEAL_II_PACKAGE_VERSION}  found at
#      ${DEAL_II_PATH}
#
#  CMAKE_BUILD_TYPE:          ${CMAKE_BUILD_TYPE}
#
#  You can now run
#      ${_make_command}                - to compile and link ${TARGET}
#      ${_make_command} debug          - to switch the build type to 'Debug'
#      ${_make_command} release        - to switch the build type to 'Release'
#      ${_make_command} clean          - to remove the generated executable as well as
#                               all intermediate compilation files
#      ${_make_command} distclean      - to clean the directory from all generated
#                               files (includes clean, runclean and the removal
#                               of the generated build system)
#      ${_make_command} info           - to view this message again
\")")

# Provide "info" target
ADD_CUSTOM_TARGET(info
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_usage.cmake
  )

# Depending on whether we link statically or allow for shared libs,
# we can or can not load plugins via external shared libs. Pass this
# down during compilation so we can disable it in the code
SET(ELASPECT_USE_SHARED_LIBS ON CACHE BOOL "If ON, we support loading shared plugin files.")
IF (DEAL_II_STATIC_EXECUTABLE STREQUAL "ON")
  MESSAGE(STATUS "Creating a statically linked executable")
  SET(ELASPECT_USE_SHARED_LIBS OFF CACHE BOOL "" FORCE)
ENDIF()

INCLUDE (CheckCXXSourceCompiles)

SET(_backup_libs ${CMAKE_REQUIRED_LIBRARIES})
LIST(APPEND CMAKE_REQUIRED_LIBRARIES ${CMAKE_DL_LIBS})
CHECK_CXX_SOURCE_COMPILES("
#include <cstddef>
#include <dlfcn.h>

int main()
{
  void *handle = dlopen (\"somelib.so\", RTLD_LAZY);
  return handle == NULL || dlerror();
}
" HAVE_DLOPEN)
SET(CMAKE_REQUIRED_LIBRARIES ${_backup_libs})

IF (NOT HAVE_DLOPEN)
  MESSAGE(STATUS "dlopen() test failed, disabling dynamic plugin loading")
  SET(ELASPECT_USE_SHARED_LIBS OFF CACHE BOOL "" FORCE)
ENDIF()

IF (ELASPECT_USE_SHARED_LIBS)
  MESSAGE(STATUS "Enabling dynamic loading of plugins from the input file")
  FOREACH(_source_file ${TARGET_SRC})
    SET_PROPERTY(SOURCE ${_source_file}
      APPEND PROPERTY COMPILE_DEFINITIONS ELASPECT_USE_SHARED_LIBS=1)
  ENDFOREACH()
ELSE()
  MESSAGE(STATUS "Disabling dynamic loading of plugins from the input file")
  FOREACH(_source_file ${TARGET_SRC})
    SET_PROPERTY(SOURCE ${_source_file}
      APPEND PROPERTY COMPILE_DEFINITIONS ELASPECT_USE_SHARED_LIBS=0)
  ENDFOREACH()
ENDIF()

# See whether we can verify that every plugin we load is compiled against
# the same deal.II library
SET(ELASPECT_HAVE_LINK_H ON CACHE BOOL "If ON, link.h exists and is usable.")
INCLUDE (CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX ("link.h" _HAVE_LINK_H)
IF (NOT _HAVE_LINK_H)
  SET(ELASPECT_HAVE_LINK_H OFF CACHE BOOL "" FORCE)
ENDIF()
IF (ELASPECT_HAVE_LINK_H)
  MESSAGE(STATUS "Enabling checking of compatible deal.II library when loading plugins")
  FOREACH(_source_file ${TARGET_SRC})
    SET_PROPERTY(SOURCE ${_source_file}
      APPEND PROPERTY COMPILE_DEFINITIONS ELASPECT_HAVE_LINK_H=1)
  ENDFOREACH()
ENDIF()


# NETCDF (c including parallel)
SET(ELASPECT_WITH_NETCDF OFF CACHE BOOL "Check if the user wants to compile elASPECT with the NETCDF libraries.")

MESSAGE(STATUS "Using ELASPECT_WITH_NETCDF = '${ELASPECT_WITH_NETCDF}'")

IF(ELASPECT_WITH_NETCDF)
  FIND_PACKAGE(NETCDF)
  IF(${NETCDF_FOUND})
    MESSAGE(STATUS "  NETCDF_INCLUDE_DIR: ${NETCDF_INCLUDE_DIR}")
    MESSAGE(STATUS "  NETCDF_LIBRARY: ${NETCDF_LIBRARY}")
    MESSAGE(STATUS :  NETCDF_VERSION: ${NETCDF_VERSION})

    INCLUDE_DIRECTORIES(${NETCDF_INCLUDE_DIRS})
    FOREACH(_T ${TARGETS})
      TARGET_LINK_LIBRARIES(${_T} ${NETCDF_LIBRARIES})
    ENDFOREACH()
  ELSE()
    MESSAGE(FATAL_ERROR "NETCDF not found. Disable ELASPECT_WITH_NETCDF or specify a hint to your installation directory with NETCDF_DIR")
  ENDIF()
ENDIF()


IF(ELASPECT_BUILD_LIBRARY)
  ADD_LIBRARY(${TARGET} ${TARGET_SRC})
  SET_TARGET_PROPERTIES(${TARGET} PROPERTIES VERSION ${ELASPECT_PACKAGE_VERSION})
ELSE()
  ADD_EXECUTABLE(${TARGET} ${TARGET_SRC})
ENDIF()
DEAL_II_SETUP_TARGET(${TARGET})

IF (ELASPECT_USE_SHARED_LIBS)
  # some systems need to explicitly link to some libraries to use dlopen
  TARGET_LINK_LIBRARIES(elaspect ${CMAKE_DL_LIBS})
ENDIF()


# Check if we can raise floating point exceptions.
SET(ELASPECT_USE_FP_EXCEPTIONS ON CACHE BOOL "If ON, floating point exception are raised in debug mode.")

IF (ELASPECT_USE_FP_EXCEPTIONS)
  INCLUDE(${CMAKE_SOURCE_DIR}/cmake/fpe_check.cmake)

  IF (HAVE_FP_EXCEPTIONS)
    MESSAGE(STATUS "Runtime floating point checks enabled.")
  ELSE()
    SET(ELASPECT_USE_FP_EXCEPTIONS OFF CACHE BOOL "" FORCE)
    MESSAGE(STATUS "No support for feenableexcept(), disabling runtime floating point exception checks.")
  ENDIF()
ENDIF()


#
## installation
#
# binary:
IF(ELASPECT_BUILD_LIBRARY)
  INSTALL(TARGET elaspect
    LIBRARY DESTINATION lib)
ELSE()
  INSTALL(TARGETS elaspect
    RUNTIME DESTINATION bin
    COMPONENT runtime)
ENDIF()

# make sure we have the rpath to our dependencies set:
SET_PROPERTY(TARGET elaspect PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

# headers:
INSTALL(DIRECTORY include/
  DESTINATION include
  COMPONENT includes
  FILES_MATCHING PATTERN "*.h")

# cmake stuff:
INSTALL(FILES ${CMAKE_BINARY_DIR}/forinstall/ElaspectConfig.cmake ${CMAKE_BINARY_DIR}/ElaspectConfigVersion.cmake
        DESTINATION "lib/cmake/Elaspect/")
