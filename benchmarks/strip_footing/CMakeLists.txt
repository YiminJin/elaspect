CMAKE_MINIMUM_REQUIRED(VERSION 3.1.0)

FIND_PACKAGE(Elaspect 1.0.0 QUIET HINTS ${Elaspect_DIR} ../ ../../ $ENV{ELASPECT_DIR})

IF (NOT Elaspect_FOUND)
  MESSAGE(FATAL_ERROR "\n"
    "Could not find a valid ELASPECT build/installation directory. "
    "Please specify the directory where you are building ELASPECT by passing\n"
    "   -D Elaspect_DIR=<path to ELASPECT>\n"
    "to cmake or by setting the environment variable ELASPECT_DIR in your shell "
    "before calling cmake. See the section 'How to write a plugin' in the "
    "manual for more information.")
ENDIF ()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

SET(TARGET "strip_footing")
PROJECT(${TARGET})

ADD_LIBRARY(${TARGET} SHARED strip_footing.cc)
ELASPECT_SETUP_PLUGIN(${TARGET})
