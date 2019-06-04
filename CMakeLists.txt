cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(project1F)
SET(VTK_DIR /Users/sam/VTK/build)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

include_directories(./include/)

add_executable(project1F src/main.cxx src/Matrix.cxx src/reader.cxx src/Camera.cxx include/Matrix.h include/reader.h include/Camera.h include/TriangleScreen.h)

if(VTK_LIBRARIES)
  target_link_libraries(project1F ${VTK_LIBRARIES})
else()
  target_link_libraries(project1F vtkHybrid)
endif()

