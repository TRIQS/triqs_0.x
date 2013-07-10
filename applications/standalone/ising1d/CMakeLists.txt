list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})
cmake_minimum_required(VERSION 2.8)
project(myproj CXX)
set(CMAKE_BUILD_TYPE Release)

find_package(TRIQS)

add_executable(ising1d ising1d.cpp)
add_definitions(-DMCTOOLS_EXPERIMENTAL)

MESSAGE ( "${TRIQS_INCLUDE_ALL}") 

include_directories(${TRIQS_INCLUDE_ALL})
target_link_libraries(ising1d ${TRIQS_LIBRARY_ALL})
