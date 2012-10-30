cmake_minimum_required(VERSION 2.8)

project(VTKHelpers)
SET(CMAKE_INCLUDE_CURRENT_DIR ON)

# For the use of shared_ptr
if(UNIX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
endif(UNIX)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

add_library(VTKHelpers VTKHelpers.cpp)
target_link_libraries(VTKHelpers ${VTK_LIBRARIES})

# Add a CMake variable so the user can choose to build the tests
option(VTKHelpers_BuildTests "VTKHelpers_BuildTests")

# Allow inheritance of include directories
if(NOT ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  set(VTKHelpers_includes ${VTKHelpers_includes} ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
  set(VTKHelpers_libraries ${VTKHelpers_libraries} VTKHelpers PARENT_SCOPE)
endif()

if(VTKHelpers_BuildTests)
  add_subdirectory(Tests)
endif(VTKHelpers_BuildTests)
