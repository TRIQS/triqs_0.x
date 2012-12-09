# Add C11 flag
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

# Compile and link with python 
# WHY ??? THIS SHOULD BE REMOVED !!!!
# Move this up : if PYTYON_SUPPRT, then EVERYONE is linked/inclu wiht ptyon...
link_libraries(${PYTHON_LIBRARY})
include_directories(${PYTHON_INCLUDE_DIRS} ${PYTHON_NUMPY_INCLUDE_DIR})

# all .cpp files...
FILE(GLOB_RECURSE ALLSOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
FOREACH(CPP ${ALLSOURCES})
 if    (NOT CPP MATCHES "/test/|/examples/|/OLD/|no-built/") 
  #MESSAGE(STATUS " lib triqs : I will compile ${CPP} ")
  set_property (GLOBAL APPEND PROPERTY TRIQS_LIB_CPP_SOURCES_LIST "${CMAKE_CURRENT_SOURCE_DIR}/${CPP}")
 endif  (NOT CPP MATCHES "/test/|/examples/|/OLD/|no-built/")
ENDFOREACH(CPP ${ALLSOURCES})
get_property(ALLSOURCES GLOBAL PROPERTY TRIQS_LIB_CPP_SOURCES_LIST)

# REMOVE THIS FROM THE CODE ...
add_definitions(-DUSE_PADE)

# MOVE THIS UP 
find_package(GMP)
SET (TRIQS_LINK_LIBS ${TRIQS_LINK_LIBS} ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})

option(Build_Triqs_General_Tools_Test "Build the test of the general tool ? ( array class , ... ) " ON)

add_all_subdirectories_with_cmakelist()

# The lib will be build from the sources
add_library(triqs ${ALLSOURCES})
target_link_libraries (triqs ${TRIQS_LINK_LIBS} )

# Install library
install(TARGETS triqs DESTINATION lib)

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


