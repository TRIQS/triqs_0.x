# boost libraries
# python and mpi have already be set before.
# here add the other boost part we want to compile

# we can not be here if boost is installed.. just check
if(NOT BOOST_SOURCE_DIR)
 message( FATAL_ERROR " cmake error in boost")
endif(NOT BOOST_SOURCE_DIR)

set(BOOST_INCLUDE_DIR ${BOOST_INSTALL_ROOT_DIR}/include  PARENT_SCOPE)

MESSAGE(STATUS "I am using boost sources from ${BOOST_SOURCE_DIR} to compile a mini boost_for_triqs ")

# what to compile ?
SET( BUILD_BOOST_SERIALIZATION ON)
#SET( BUILD_BOOST_FILESYSTEM ON)
SET( BUILD_BOOST_SYSTEM ON)

SET( BUILD_BOOST_FILESYSTEM OFF)
SET (BUILD_BOOST_THREAD OFF)
SET (BUILD_BOOST_DATE_TIME OFF )
SET (BUILD_BOOST_PROGRAM_OPTIONS OFF)
SET (BUILD_BOOST_REGEX OFF) 

set(ALPS_BOOST_LIBRARY_NAME "boost_for_triqs" )
include (CompileBoost)

include_directories(${BOOST_SOURCE_DIR}) # this one is for the compilation in this dir
if(TRIQS_BUILD_STATIC)
  set(BOOST_LIBRARY boost_for_triqs mpi ${ALPS_EXTRA_LIBRARIES} PARENT_SCOPE )
else(TRIQS_BUILD_STATIC)
  set(BOOST_LIBRARY boost_for_triqs ${ALPS_EXTRA_LIBRARIES} PARENT_SCOPE )
endif(TRIQS_BUILD_STATIC)
set(BOOST_INCLUDE_DIR ${BOOST_SOURCE_DIR} PARENT_SCOPE) # this is for later compilation in the project

# install the python sources 
SET(PYTHON_SOURCES 
 ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py 
 )

if(TRIQS_BUILD_STATIC)

  set_property(GLOBAL APPEND PROPERTY Modules_actually_compiled mpi)

else(TRIQS_BUILD_STATIC)

 install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/boost)
 install (FILES ${CMAKE_CURRENT_BINARY_DIR}/mpi.so  DESTINATION ${TRIQS_PYTHON_LIB_DEST}/boost)

 # for tests
 FILE  (COPY ${PYTHON_SOURCES} DESTINATION ${CMAKE_BINARY_DIR}/pytriqs/boost)
 execute_process(COMMAND ln -fs ${CMAKE_CURRENT_BINARY_DIR}/mpi.so  ${CMAKE_BINARY_DIR}/pytriqs/boost)

endif(TRIQS_BUILD_STATIC)
