
option(BUILD_BCPLUGIN              "Building BCPlugin" OFF)
option(BUILD_BCPLUGIN_WITH_SICONOS "Building BCPlugin with Siconos" OFF)

if(NOT BUILD_BCPLUGIN)
  return()
endif()

include(cmake/Findsiconos.cmake)

message (STATUS "USE_SICONOS=${BUILD_BCPLUGIN_WITH_SICONOS}")

set(target CnoidBCPlugin)

set(sources
  BCPlugin.cpp
  BCSimulatorItem.cpp
  BCConstraintForceSolver.cpp
  BCCoreSiconos.cpp
  BCCoreQMR.cpp
  )

set(headers
  BCSimulatorItem.h
  BCConstraintForceSolver.h
  BCCoreSiconos.h
  BCCoreQMR.h
  )

if(BUILD_BCPLUGIN_WITH_SICONOS)
  include_directories(
    ${SICONOS_INCLUDE_DIR}
  )
  add_definitions( -DBUILD_BCPLUGIN_WITH_SICONOS )
endif()

make_gettext_mofiles(${target} mofiles)
add_cnoid_plugin(${target} SHARED ${sources} ${headers} ${mofiles})
target_link_libraries(${target} CnoidBodyPlugin)
if(BUILD_BCPLUGIN_WITH_SICONOS)
  target_link_libraries(${target} siconos_numerics)
endif()
apply_common_setting_for_plugin(${target} "${headers}")

if(ENABLE_PYTHON)
#  add_subdirectory(python)
endif()
