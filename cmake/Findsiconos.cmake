# - Try to find siconos
# Once done this will define
#
#  SICONOS_FOUND - system has siconos
#  SICONOS_INCLUDE_DIR - the siconos include directory
#  SICONOS_LIBRARY - Link these to use siconos
#

FIND_PATH(SICONOS_INCLUDE_DIR NAMES SiconosKernel.hpp
  PATHS
  $ENV{SICONOS_INCLUDE_DIR}
  /usr/include/siconos
  /usr/local/include/siconos
)

FIND_LIBRARY(SICONOS_LIBRARY NAMES siconos_kernel
  PATHS
  $ENV{SICONOS_LIB_DIR}
  /usr/lib
  /usr/local/lib
)

IF(SICONOS_INCLUDE_DIR AND SICONOS_LIBRARY)
   SET(SICONOS_FOUND TRUE)
ENDIF(SICONOS_INCLUDE_DIR AND SICONOS_LIBRARY)

# show the SICONOS_INCLUDE_DIR and SICONOS_LIBRARY variables only in the advanced view
IF(SICONOS_FOUND)
  MARK_AS_ADVANCED(SICONOS_INCLUDE_DIR SICONOS_LIBRARY )
ENDIF(SICONOS_FOUND)
