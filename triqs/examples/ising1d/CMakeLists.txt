cmake_minimum_required(VERSION 2.8)

project(myproj CXX)

set(CMAKE_BUILD_TYPE Release)

set(TRIQS_PATH /path/to/installed/triqs/base)
include(${TRIQS_PATH}/share/triqs/cmake/TRIQSConfig.cmake)

add_executable(ising1d ising1d.cpp)
add_definitions(-DMCTOOLS_EXPERIMENTAL)

include_directories(${TRIQS_INCLUDE} ${EXTRA_INCLUDE} ${CBLAS_INCLUDE} ${FFTW_INCLUDE})
target_link_libraries(ising1d ${TRIQS_LIBRARY})
