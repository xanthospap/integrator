cmake_minimum_required(VERSION 3.8)

project(
  orbit_integrator
  VERSION 1.1.0
  DESCRIPTION "Satellite Orbit Integrators for POD"
  LANGUAGES CXX
)

# Ensure required libraries are available
find_package(datetime REQUIRED)
find_package(iers REQUIRED)

# Pass the library dependencies to subdirectories
set(PROJECT_DEPENDENCIES iers datetime)

# Define an option for building binaries (defaults to ON)
option(BUILD_PROGS "Enable building of binaries" ON)
if(BUILD_PROGS)
  find_package(sp3 REQUIRED)
endif()

# compiler flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED On)
set(CMAKE_CXX_EXTENSIONS Off)
add_compile_options(-Wall
  -Wextra
  -Werror
  -pedantic
  -W
  -Wshadow
  $<$<CONFIG:Release>:-O2>
  $<$<CONFIG:Release>:-march=native>
  $<$<CONFIG:Debug>:-g>
  $<$<CONFIG:Debug>:-pg>
  $<$<CONFIG:Debug>:-Wdisabled-optimization>
  )
add_compile_definitions(
  $<$<CONFIG:Debug>:DEBUG>
)

# the library and includes
add_library(integrator)
target_include_directories(integrator
PUBLIC
  # Public headers during build time
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
  # Public headers during install time
  $<INSTALL_INTERFACE:include/orbit_integrator>
PRIVATE
  # Private/internal headers only during build time (src/core)
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  # Internal/private headers during install time (from 
  # include/orbit_integrator/core)
  # $<INSTALL_INTERFACE:include/orbit_integrator/core>
)

# library source code
add_subdirectory(src)

if(BUILD_PROGS)
  add_subdirectory(bin)
endif()

# install library
install(TARGETS integrator
        EXPORT integratorTargets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin
)

# Install public headers at: $PREFIX/orbit_integrator/...
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/
	DESTINATION include/orbit_integrator
)

# Install private headers to a separate directory $PREFIX/datetime/core
# install(DIRECTORY ${CMAKE_SOURCE_DIR}/src/core/
#         DESTINATION include/datetime/core
# )

install(EXPORT integratorTargets
        FILE integratorTargets.cmake
        NAMESPACE dso::
        DESTINATION lib/cmake
)

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  "integratorConfigVersion.cmake"
  VERSION ${integrator_version}
  COMPATIBILITY AnyNewerVersion
)

install(FILES "integratorConfig.cmake" "${CMAKE_CURRENT_BINARY_DIR}/integratorConfigVersion.cmake"
  DESTINATION lib/cmake
)
