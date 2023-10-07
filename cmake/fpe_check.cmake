# Copyright (C) 2023 by Yimin Jin.
#
#  This file is part of elASPECT.
#
#  elASPECT is modified from the free software ASPECT; you can 
#  redistribute it and/or modify it under the terms of the GNU 
#  General Public License as published by the Free Software 
#  Foundation; either version 2, or (at your option) any later 
#  version.
#
#  elASPECT is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with elASPECT; see the file LICENSE.  If not see
#  <http://www.gnu.org/licenses/>.

INCLUDE (CheckCXXSourceRuns)

SET(_backup_flags ${CMAKE_REQUIRED_FLAGS})
SET(_backup_libs ${CMAKE_REQUIRED_LIBRARIES})
SET(_backup_includes ${CMAKE_REQUIRED_INCLUDES})

SET(_build "RELEASE")
STRING(TOLOWER "${CMAKE_BUILD_TYPE}" _cmake_build_type)
IF("${_cmake_build_type}" MATCHES "debug")
  SET(_build "DEBUG")
ENDIF()

LIST(APPEND CMAKE_REQUIRED_FLAGS "${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${_build}}")
LIST(APPEND CMAKE_REQUIRED_LIBRARIES ${DEAL_II_TARGET_${_build}})
LIST(APPEND CMAKE_REQUIRED_INCLUDES ${DEAL_II_INCLUDE_DIRS})

CHECK_CXX_SOURCE_RUNS("
#include <fenv.h>
#include <limits>
#include <sstream>

#include <deal.II/base/utilities.h>

int main()
{
  // Some implementations seem to not initialize the FPE bits to zero.
  // Make sure we start from a clean state
  feclearexcept(FE_DIVBYZERO|FE_INVALID);

  // Enable floating point exceptions
  feenableexcept(FE_DIVBYZERO|FE_INVALID);

  std::ostringstream description;
  const double lower_bound = -std::numeric_limits<double>::max();

  description << lower_bound;
  description << dealii::Utilities::string_to_int (\"1\");

  return 0;
}
" HAVE_FP_EXCEPTIONS)

SET(CMAKE_REQUIRED_FLAGS ${_backup_flags})
SET(CMAKE_REQUIRED_LIBRARIES ${_backup_libs})
SET(CMAKE_REQUIRED_INCLUDES ${_backup_includes})
