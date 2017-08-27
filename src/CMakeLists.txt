project (aither)

# set version number
set (aither_VERSION_MAJOR 0)
set (aither_VERSION_MINOR 8)
set (aither_VERSION_PATCH 0)

message (STATUS "Cmake version ${CMAKE_VERSION}")
message (STATUS "Configuring aither v${aither_VERSION_MAJOR}.${aither_VERSION_MINOR}.${aither_VERSION_PATCH}")
message (STATUS "Installation directory: ${CMAKE_INSTALL_PREFIX}")
message (STATUS "Include path: ${CMAKE_INCLUDE_PATH}")
message (STATUS "Library path: ${CMAKE_LIBRARY_PATH}")
message (STATUS "Prefix path: ${CMAKE_PREFIX_PATH}")
message (STATUS "Module path: ${CMAKE_MODULE_PATH}")
message (STATUS "Build type is: ${CMAKE_BUILD_TYPE}")
message (STATUS "System name: ${CMAKE_SYSTEM_NAME}")

# add include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${PROJECT_BINARY_DIR})

# configure header file to pass cmake variables
configure_file (
  "${CMAKE_SOURCE_DIR}/include/macros.hpp.in"
  "${PROJECT_BINARY_DIR}/macros.hpp"
  )

# set all source files
set(sources
  main.cpp
  boundaryConditions.cpp
  eos.cpp
  fluid.cpp
  fluxJacobian.cpp
  genArray.cpp
  input.cpp
  inputStates.cpp
  inviscidFlux.cpp
  kdtree.cpp
  matrix.cpp
  output.cpp
  parallel.cpp
  plot3d.cpp
  primVars.cpp
  procBlock.cpp
  range.cpp
  resid.cpp
  slices.cpp
  source.cpp
  thermodynamic.cpp
  transport.cpp
  turbulence.cpp
  uncoupledScalar.cpp
  utility.cpp
  viscousFlux.cpp
  wallData.cpp
  wallLaw.cpp
  )

# create targets for executable and libraries
add_executable (aither ${sources})
add_library (aitherShared SHARED ${sources})
add_library (aitherStatic STATIC ${sources})

# rename libraries
set_target_properties (aitherStatic PROPERTIES OUTPUT_NAME aither)
set_target_properties (aitherShared PROPERTIES OUTPUT_NAME aither)

# set c++ standard for all targets
set_property (TARGET aither aitherStatic aitherShared PROPERTY CXX_STANDARD 14)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
message (STATUS "Using C++14")

# add compiler specific flags
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  message(STATUS "Changing stdlib for Clang")
  set(COMPILER_SPEC_FLAGS "-stdlib=libc++")
endif()

# set additional c++ flags for all build types
if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMPILER_SPEC_FLAGS} /O2 /DNDEBUG")
else()
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMPILER_SPEC_FLAGS} -Wall -pedantic -march=native -O3 -DNDEBUG")
endif()

# release build type
if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
   set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${COMPILER_SPEC_FLAGS}")
else()
  set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${COMPILER_SPEC_FLAGS} -Wall -pedantic -march=native")
endif()

# debug build type
if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
   set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${COMPILER_SPEC_FLAGS} /Wall /Od /PROFILE")
else()
   set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${COMPILER_SPEC_FLAGS} -Wall -pedantic -O0 -ggdb -pg --coverage")
endif()

# create profile build type
if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
   set (CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS} ${COMPILER_SPEC_FLAGS} /PROFILE")
else()
   set (CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS} ${COMPILER_SPEC_FLAGS} -pg")
endif()

# set compiler flags based on release type
if (CMAKE_BUILD_TYPE MATCHES "release")
   message (STATUS "C++ flags: " ${CMAKE_CXX_FLAGS_RELEASE})
elseif (CMAKE_BUILD_TYPE MATCHES "debug")
   message (STATUS "C++ flags: " ${CMAKE_CXX_FLAGS_DEBUG})
elseif (CMAKE_BUILD_TYPE MATCHES "profile")
   message (STATUS "C++ flags: " ${CMAKE_CXX_FLAGS_PROFILE})
elseif (CMAKE_BUILD_TYPE MATCHES "relwithdebinfo")
   message (STATUS "C++ flags: " ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
elseif (CMAKE_BUILD_TYPE MATCHES "minsizerel")
   message (STATUS "C++ flags: " ${CMAKE_CXX_FLAGS_MINSIZEREL})
else ()
   message (STATUS "C++ flags: " ${CMAKE_CXX_FLAGS})
endif ()

# determine if MPI_DIR has been set
if (MPI_DIR)
   message (STATUS "Looking for MPI here: ${MPI_DIR}")
   set (CMAKE_PREFIX_PATH ${MPI_DIR})
endif ()

# get mpi and link to all targets
find_package (MPI REQUIRED)
include_directories (SYSTEM "${MPI_INCLUDE_PATH}")
target_link_libraries (aither ${MPI_C_LIBRARIES})
target_link_libraries (aitherStatic ${MPI_C_LIBRARIES})
target_link_libraries (aitherShared ${MPI_C_LIBRARIES})

# install executable, libraries, and includes
install (TARGETS aither aitherStatic aitherShared
	ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib
	RUNTIME DESTINATION bin)
FILE (GLOB files "${CMAKE_SOURCE_DIR}/include/*.hpp")
install (FILES ${files} DESTINATION include)
install (FILES "${PROJECT_BINARY_DIR}/macros.hpp" DESTINATION include)
