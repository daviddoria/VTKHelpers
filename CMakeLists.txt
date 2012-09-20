cmake_minimum_required(VERSION 2.8)

project(VTKHelpers)
SET(CMAKE_INCLUDE_CURRENT_DIR ON)

# For the use of shared_ptr
if(UNIX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif(UNIX)

# VTK
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})
set(VTKHelpers_libraries ${VTKHelpers_libraries} ${VTK_LIBRARIES})

# Give the compiler all of the required include directories
include_directories(${VTKHelpers_include_dirs})

# Create the library
add_library(VTKHelpers VTKHelpers.cpp)
target_link_libraries(VTKHelpers ${VTK_LIBRARIES})
set(VTKHelpers_libraries ${VTKHelpers_libraries} VTKHelpers) # Add the library to the list of libraries used

# Add a CMake variable so the user can choose to build the tests
SET(VTKHelpers_BuildTests ON CACHE BOOL "VTKHelpers_BuildTests")

# Allow this project to be detected and used as a submodule
set(VTKHelpers_include_dirs ${VTKHelpers_include_dirs} ${CMAKE_CURRENT_SOURCE_DIR})
set_property(GLOBAL PROPERTY VTKHelpersIncludeDirs ${VTKHelpers_include_dirs})
set_property(GLOBAL PROPERTY VTKHelpersLibs ${VTKHelpers_libraries})

# Build the tests if requested
if(VTKHelpers_BuildTests)
  add_subdirectory(Tests)
endif(VTKHelpers_BuildTests)
