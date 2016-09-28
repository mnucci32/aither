cmake_minimum_required (VERSION 2.6)
project (aither)

# set version number
set (aither_VERSION_MAJOR 0)
set (aither_VERSION_MINOR 3)
set (aither_VERSION_PATCH 1))

include_directories(${CMAKE_SOURCE_DIR}/include)

# configure header file to pass cmake variables
configure_file (
	       "${PROJECT_SOURCE_DIR}/macros.hpp.in"
	       "${PROJECT_BINARY_DIR}/macros.hpp"
)

add_executable(aither main.cpp)

# set standard
set_property (TARGET aither PROPERTY CXX_STANDARD 14)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

# get mpi
find_package (MPI REQUIRED)
include_directories (SYSTEM "${MPI_INCLUDE_PATH}")
target_link_libraries (aither ${MPI_C_LIBRARIES})

# add install rule
install (TARGETS aither DESTINATION bin)
install (FILES "${PROJECT_BINARY_DIRECTORY}/macro.hpp" DESTINATION include)
