cmake_minimum_required(VERSION 3.10)

PROJECT(proj6)
SET(VTK_DIR /Users/hank/Hartree/VTK/install/lib/cmake/vtk-8.1)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(proj6 MACOSX_BUNDLE proj6)

SET(CMAKE_CXX_LINK_FLAGS "-framework OpenGL")
if(VTK_LIBRARIES)
  target_link_libraries(proj6 ${VTK_LIBRARIES})
else()
  target_link_libraries(proj6 vtkHybrid)
endif()


