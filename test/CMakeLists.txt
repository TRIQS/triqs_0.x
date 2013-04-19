# some specific flags for tests...
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
include_directories( ${CMAKE_SOURCE_DIR} )

SET( link_libs ${LAPACK_LIBS}  ${BOOST_LIBRARY} ${ALPS_EXTRA_LIBRARIES})
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
 list (REMOVE_DUPLICATES link_libs)
ENDIF( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
link_libraries( ${link_libs} triqs_utility )

if (CMAKE_COMPILER_IS_CLANG) 
 #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wunused-variable")
endif()

add_subdirectory(triqs)
if (PythonSupport)
 add_subdirectory(pytriqs)
endif()
add_subdirectory(speed)

