# To allow includes like "VTKHelpers.h" instead of "VTKHelpers/VTKHelpers.h"
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

add_executable(TestVTKHelpers TestVTKHelpers.cpp)
target_link_libraries(TestVTKHelpers VTKHelpers)
