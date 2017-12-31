cmake_minimum_required(VERSION 2.8)
set( CMAKE_VERBOSE_MAKEFILE on )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )

project(quiz-2d-KF-Laser)
include_directories( ./ ./include )
add_executable(${PROJECT_NAME} src/kalman_filter.cpp  src/main.cpp  src/tracking.cpp)

