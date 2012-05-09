cmake_minimum_required(VERSION 2.8)

project(VTKHelpers)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# This is needed for shared_ptr and the trick using enable_if and if_fundamental to allow scalars
# to be treated as the 0th component of a vector.
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x")

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

add_library(libVTKHelpers VTKHelpers.cpp)
target_link_libraries(libVTKHelpers ${VTK_LIBRARIES})

if(BuildTests)
  add_subdirectory(Tests)
endif(BuildTests)

