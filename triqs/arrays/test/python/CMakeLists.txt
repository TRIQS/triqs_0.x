# Name of extended python interpreter to be compiled
SET ( TestList pythonArrayTest )
 
# Tool to help to include the wrappers in main and declare the built-in module
#Find_package(PythonModuleEmbeddingTool)

# Compile and link with python 
#include_directories(${PYTHON_INCLUDE_DIRS} ${PYTHON_NUMPY_INCLUDE_DIR})

# We use the headers of triqs and Base in submodules + Blitz + FFTW
# ${CMAKE_CURRENT_BINARY_DIR} is for the main_includewrap.cpp main_appendinittab.cpp files
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Compile main and link the whole thing....
SET(PythonPathAdd " ${CMAKE_INSTALL_PREFIX}/lib/pytriqs  ")
add_definitions( 
  -DPYTHONPATH_ADD=${PythonPathAdd}
  -DTRIQS_ARRAYS_WITH_PYTHON_SUPPORT
  -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  )

# Add the necessary libs to the main target.
link_libraries( 
 ${PYTHON_LIBRARY} ${PYTHON_EXTRA_LIBS} 
  #  ${Boost_LIBRARIES} 
  ${LAPACK_LIBS} 
 ${BOOST_LIBRARY} ${ALPS_EXTRA_LIBRARIES}
)


find_package(TestCSimple)

FOREACH( TestName  ${TestList} )
 add_executable( ${TestName}  ${CMAKE_CURRENT_SOURCE_DIR}/${TestName}.cpp )
 add_test_C_simple( ${TestName} )
 if (MPI_FOUND)
  if( BOOST_SOURCE_DIR)
   if (TRIQS_BUILD_STATIC)
    target_link_libraries( ${TestName}  ${CMAKE_BINARY_DIR}/boost/mpi.a)
   else (TRIQS_BUILD_STATIC)
    #target_link_libraries( ${TestName}  ${CMAKE_BINARY_DIR}/boost/mpi.so)
   endif (TRIQS_BUILD_STATIC)
  endif( BOOST_SOURCE_DIR)
 endif(MPI_FOUND)
ENDFOREACH( TestName  ${TestList} )

