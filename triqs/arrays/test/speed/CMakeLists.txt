FILE(GLOB TestList RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

enable_testing()

include_directories( ${CMAKE_CURRENT_BINARY_DIR})

link_libraries( 
 ${Boost_LIBRARIES} 
 ${LAPACK_LIBS} 
 # ${HDF5_LIBRARIES}
 ${ALPS_LIBRARIES}   
 
 )

link_libraries( ${link_libs} triqs_utility ) # adding triqs for the stack_trace fun used in exceptions... 

find_package(TestCSimple)

FOREACH( TestName1  ${TestList} )
 STRING(REPLACE ".cpp" "" TestName ${TestName1})
 add_executable( ${TestName}  ${CMAKE_CURRENT_SOURCE_DIR}/${TestName} )
 #add_test_C_simple( ${TestName} )
ENDFOREACH( TestName1  ${TestList} )



