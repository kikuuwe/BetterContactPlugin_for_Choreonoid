
option(BUILD_BCPLUGIN              "Building BCPlugin" OFF)
option(BUILD_BCPLUGIN_WITH_SICONOS "Building BCPlugin with Siconos" OFF)

if(NOT BUILD_BCPLUGIN)
  return()
endif()



configure_file (${CMAKE_CURRENT_SOURCE_DIR}/Config.hpp.in
                ${CMAKE_CURRENT_BINARY_DIR}/Config.hpp)

message (STATUS "Generating Config.hpp ...")
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
  )

if(BUILD_BCPLUGIN_WITH_SICONOS)
  include_directories(
    /usr/local/include/Siconos/Numerics
  )
endif()

link_directories(/usr/local/lib/x86_64-linux-gnu)

 
make_gettext_mofiles(${target} mofiles)
add_cnoid_plugin(${target} SHARED ${sources} ${headers} ${mofiles})
target_link_libraries(${target} CnoidBodyPlugin)
apply_common_setting_for_plugin(${target} "${headers}")


if(BUILD_BCPLUGIN_WITH_SICONOS)
  target_link_libraries(${target} SiconosNumerics)
endif()


if(ENABLE_PYTHON)
#  add_subdirectory(python)
endif()
