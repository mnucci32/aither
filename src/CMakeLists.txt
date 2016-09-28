cmake_minimum_required (VERSION 2.6)
project (aither)

# set version number
set (aither_VERSION_MAJOR 0)
set (aither_VERSION_MINOR 3)
set (aither_VERSION_PATCH 1)


message(STATUS "CMAKE_SOURCE_DIR" ${CMAKE_SOURCE_DIR})
message(STATUS "CMAKE_CURRENT_SOURCE_DIR" ${CMAKE_CURRENT_SOURCE_DIR})
message(STATUS "PROJECT_SOURCE_DIR" ${PROJECT_SOURCE_DIR})
message(STATUS "PROJECT_BINARY_DIR" ${PROJECT_BINARY_DIR})

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${PROJECT_BINARY_DIR})

# configure header file to pass cmake variables
configure_file (
  "${CMAKE_SOURCE_DIR}/include/macros.hpp.in"
  "${PROJECT_BINARY_DIR}/macros.hpp"
  )

add_executable(aither
  main.cpp
  boundaryConditions.cpp
  eos.cpp
  fluxJacobian.cpp
  genArray.cpp
  input.cpp
  inviscidFlux.cpp
  kdtree.cpp
  matrix.cpp
  output.cpp
  parallel.cpp
  plot3d.cpp
  primVars.cpp
  procBlock.cpp
  resid.cpp
  slices.cpp
  source.cpp
  turbulence.cpp
  uncoupledScalar.cpp
  utility.cpp
  viscousFlux.cpp
  )

# set standard
set_property (TARGET aither PROPERTY CXX_STANDARD 14)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

# set compiler options
add_compile_options(-Wall -pedantic -march=native -O3)

# get mpi
find_package (MPI REQUIRED)
include_directories (SYSTEM "${MPI_INCLUDE_PATH}")
target_link_libraries (aither ${MPI_C_LIBRARIES})

# add install rule
install (TARGETS aither DESTINATION bin)
FILE(GLOB files "${CMAKE_SOURCE_DIR}/include/*.hpp")
install (FILES ${files} DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/macros.hpp" DESTINATION include)
