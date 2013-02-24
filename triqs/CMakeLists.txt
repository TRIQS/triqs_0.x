# Add C11 flag
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

# Compile and link with python 
# WHY ??? THIS SHOULD BE REMOVED !!!!
# Move this up : if PYTYON_SUPPRT, then EVERYONE is linked/inclu wiht ptyon...
link_libraries(${PYTHON_LIBRARY})
include_directories(${PYTHON_INCLUDE_DIRS} ${PYTHON_NUMPY_INCLUDE_DIR})

# all .cpp files...
FILE(GLOB_RECURSE ALLSOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
SET( EXCLUSION_LIST "/test/|/examples/|/OLD/|no-built/")
FOREACH(CPP ${ALLSOURCES})
 if    (NOT CPP MATCHES ${EXCLUSION_LIST})
 if  (CPP MATCHES "utility|mc_tools") 
  set_property (GLOBAL APPEND PROPERTY TRIQS_UTIL_LIB_CPP_SOURCES_LIST "${CMAKE_CURRENT_SOURCE_DIR}/${CPP}")
 else  (CPP MATCHES "utility|mc_tools")
  #MESSAGE(STATUS " lib triqs : I will compile ${CPP} ")
  set_property (GLOBAL APPEND PROPERTY TRIQS_LIB_CPP_SOURCES_LIST "${CMAKE_CURRENT_SOURCE_DIR}/${CPP}")
 endif (CPP MATCHES "utility|mc_tools")
 endif  (NOT CPP MATCHES ${EXCLUSION_LIST})
ENDFOREACH(CPP ${ALLSOURCES})
get_property(ALLSOURCES GLOBAL PROPERTY TRIQS_LIB_CPP_SOURCES_LIST)
get_property(ALLSOURCES_UTIL GLOBAL PROPERTY TRIQS_UTIL_LIB_CPP_SOURCES_LIST)

# MOVE THIS UP 
find_package(GMP)
SET (TRIQS_LINK_LIBS ${TRIQS_LINK_LIBS} ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})

option(Build_Triqs_General_Tools_Test "Build the test of the general tool ? ( array class , ... ) " ON)

add_all_subdirectories_with_cmakelist()

# The lib will be build from the sources
add_library(triqs_utility ${ALLSOURCES_UTIL})
target_link_libraries (triqs_utility -lrt ${TRIQS_LINK_LIBS} )

add_library(triqs ${ALLSOURCES})
target_link_libraries (triqs triqs_utility ${TRIQS_LINK_LIBS}  )

# Install library
install(TARGETS triqs DESTINATION lib)
install(TARGETS triqs_utility DESTINATION lib)

# TRIQS cmake file
configure_file(../cmake/TRIQSConfig.cmake.in Config/TRIQSConfig.cmake)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/Config/TRIQSConfig.cmake DESTINATION share/triqs/cmake)

# all .hpp files...
FILE(GLOB_RECURSE ALLHEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp)

# Install headers in good include dir
FOREACH(HEADER ${ALLHEADERS})
 get_filename_component(P ${HEADER} PATH)
 install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/${HEADER}" DESTINATION "include/triqs/${P}")
 #MESSAGE(STATUS "I will install ${HEADER} of ${P}")
ENDFOREACH(HEADER ${ALLHEADERS})


