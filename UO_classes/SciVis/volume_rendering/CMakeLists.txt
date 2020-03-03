cmake_minimum_required(VERSION 3.10)

PROJECT(volume_rendering)
SET(VTK_DIR /Users/hank/Hartree/VTK/install/lib/cmake/vtk-8.1)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(volume_rendering MACOSX_BUNDLE volume_rendering)

SET(CMAKE_CXX_LINK_FLAGS "-framework OpenGL")
if(VTK_LIBRARIES)
  target_link_libraries(volume_rendering ${VTK_LIBRARIES})
else()
  target_link_libraries(volume_rendering vtkHybrid)
endif()


