cmake_minimum_required(VERSION 2.8)

PROJECT(proj7a)
SET(VTK_DIR /Users/jweiler/VTK-build)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(proj7a proj7a)

SET(CMAKE_CXX_LINK_FLAGS "-framework OpenGL")

if(VTK_LIBRARIES)
  target_link_libraries(proj7a ${VTK_LIBRARIES})
else()
  target_link_libraries(proj7a vtkHybrid)
endif()
