FILE(GLOB TestList RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
#MESSAGE(STATUS " Tests for arrays libraries :  ${TestList}")

enable_testing()

include_directories( ${CMAKE_CURRENT_BINARY_DIR})

find_package(TestCSimple)

STRING (REPLACE "/" ";" pa ${CMAKE_CURRENT_SOURCE_DIR} )
LIST(GET pa -1 dir)
STRING(COMPARE EQUAL ${dir} "C++_with_py_support" python_support)

# First the version for pure C, then with python support
if (NOT python_support)

 add_definitions( -UTRIQS_WITH_PYTHON_SUPPORT)

SET( link_libs ${LAPACK_LIBS}  ${BOOST_LIBRARY} ${ALPS_EXTRA_LIBRARIES})
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
 list (REMOVE_DUPLICATES link_libs)
ENDIF( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

link_libraries( ${link_libs} triqs ) # adding triqs for the stack_trace fun used in exceptions... 

FOREACH( TestName1  ${TestList} )
 STRING(REPLACE ".cpp" "" TestName ${TestName1})
 add_executable( ${TestName}_nopy  ${CMAKE_CURRENT_SOURCE_DIR}/${TestName}.cpp )
 add_test_C_simple( ${TestName} _nopy )
ENDFOREACH( TestName1  ${TestList} )

else (NOT python_support)

# We use the headers of triqs and base in submodules + Blitz + FFTW
# ${CMAKE_CURRENT_BINARY_DIR} is for the main_includewrap.cpp main_appendinittab.cpp files
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Compile main and link the whole thing....
SET(PythonPathAdd " ${CMAKE_INSTALL_PREFIX}/lib/pytriqs  ")
add_definitions( 
  -DPYTHONPATH_ADD=${PythonPathAdd}
  -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  )

# Add the necessary libs to the main target.
SET( link_libs
  triqs 
 ${PYTHON_LIBRARY} ${PYTHON_EXTRA_LIBS} 
  #  ${Boost_LIBRARIES} 
  ${LAPACK_LIBS} 
 ${BOOST_LIBRARY} ${ALPS_EXTRA_LIBRARIES}
)

# does not work because of optimized
#MESSAGE( "link_libs = ${link_libs}")

#IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
# list (REMOVE_DUPLICATES link_libs)
#ENDIF( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")


link_libraries( ${link_libs}) 

FOREACH( TestName1  ${TestList} )
 STRING(REPLACE ".cpp" "" TestName ${TestName1})
 add_executable( ${TestName}_wpy  ${CMAKE_CURRENT_SOURCE_DIR}/${TestName}.cpp )
 add_test_C_simple( ${TestName} _wpy )
 if (MPI_FOUND)
  if( BOOST_SOURCE_DIR)
   if (TRIQS_BUILD_STATIC)
    target_link_libraries( ${TestName}_wpy  ${CMAKE_BINARY_DIR}/foreignlibs/boost/mpi.a)
   else (TRIQS_BUILD_STATIC)
    # why ? does not work on os X anyway
    #target_link_libraries( ${TestName}_wpy  ${CMAKE_BINARY_DIR}/boost/mpi.so)
   endif (TRIQS_BUILD_STATIC)
  endif( BOOST_SOURCE_DIR)
 endif(MPI_FOUND)
ENDFOREACH( TestName1  ${TestList} )

endif (NOT python_support)


