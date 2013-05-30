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
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
SET( LINK_RT "")
else()
SET( LINK_RT "-lrt")
endif()
target_link_libraries (triqs_utility ${LINK_RT} ${TRIQS_LINK_LIBS} )
add_library(triqs ${ALLSOURCES})
target_link_libraries (triqs triqs_utility ${TRIQS_LINK_LIBS}  )

# Install library
install(TARGETS triqs DESTINATION lib)
install(TARGETS triqs_utility DESTINATION lib)

# TRIQS cmake file
get_property(ALL_DEFS GLOBAL PROPERTY TRIQS_COMPIL_DEFINITIONS)
IF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
 configure_file(../cmake/TRIQSConfig.cmake.in Config/TRIQSConfig.cmake @ONLY)
else()
 configure_file(../cmake/TRIQSConfig.cmake.in Config/TRIQSConfig1.cmake @ONLY )
 EXEC_PYTHON_SCRIPT( "import re;open('${CMAKE_CURRENT_BINARY_DIR}/Config/TRIQSConfig.cmake','w').write(re.sub(r';*\\s*-framework\\s+([a-zA-Z_][a-zA-Z_0-9]*)',r' \"-framework \\1\"',open('${CMAKE_CURRENT_BINARY_DIR}/Config/TRIQSConfig1.cmake').read()))", ConfigFile)
endif()

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/Config/TRIQSConfig.cmake DESTINATION share/triqs/cmake)

# all .hpp files...
FILE(GLOB_RECURSE ALLHEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp)

# Install headers in good include dir
FOREACH(HEADER ${ALLHEADERS})
 get_filename_component(P ${HEADER} PATH)
 install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/${HEADER}" DESTINATION "include/triqs/${P}")
 #MESSAGE(STATUS "I will install ${HEADER} of ${P}")
ENDFOREACH(HEADER ${ALLHEADERS})


