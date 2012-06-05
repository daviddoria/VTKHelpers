cmake_minimum_required(VERSION 2.8)

project(VTKHelpers)

# For the use of shared_ptr
if(UNIX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif(UNIX)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

add_library(VTKHelpers VTKHelpers.cpp)
target_link_libraries(VTKHelpers ${VTK_LIBRARIES})

SET(VTKHelpers_BuildTests ON CACHE BOOL "VTKHelpers_BuildTests" FORCE)
if(VTKHelpers_BuildTests)
  add_subdirectory(Tests)
endif(VTKHelpers_BuildTests)

