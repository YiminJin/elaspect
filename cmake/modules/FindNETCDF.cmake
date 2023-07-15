SET(NETCDF_DIR "" CACHE PATH "An optional hint to a NETCDF installation")
SET_IF_EMPTY(NETCDF_DIR "$ENV{NETCDF_DIR}")

FIND_PATH(NETCDF_INCLUDE_DIR
  NAMES netcdf_meta.h
  HINTS ${NETCDF_DIR}
  PATH_SUFFIXES netcdf include
  )

FIND_LIBRARY(NETCDF_LIBRARY
  NAMES netcdf
  HINTS ${NETCDF_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )


SET(_header "${NETCDF_INCLUDE_DIR}/netcdf_meta.h")

SET(NETCDF_VERSION "unknown")

IF(EXISTS ${_header})
  FILE(STRINGS "${_header}" _version_line
    REGEX "#define.*NC_VERSION"
    )
  STRING(REGEX REPLACE ".*\"(.+)\"" "\\1" NETCDF_VERSION
    "${_version_line}"
    )
ENDIF()


IF(NETCDF_INCLUDE_DIR AND NETCDF_LIBRARY)
  SET(NETCDF_FOUND TRUE)
  SET(NETCDF_LIBRARIES ${NETCDF_LIBRARY})
  SET(NETCDF_INCLUDE_DIRS ${NETCDF_INCLUDE_DIR})
ELSE()
    SET(NETCDF_FOUND FALSE)
ENDIF()
