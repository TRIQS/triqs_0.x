SET ( TestList speed speed2 gemv )

enable_testing()

include_directories( ${CMAKE_CURRENT_BINARY_DIR})

link_libraries( 
 ${Boost_LIBRARIES} 
 ${LAPACK_LIBS} 
 # ${HDF5_LIBRARIES}
 ${ALPS_LIBRARIES}   
 
 )

find_package(TestCSimple)

FOREACH( TestName  ${TestList} )
 add_executable( ${TestName}  ${CMAKE_CURRENT_SOURCE_DIR}/${TestName}.cpp )
 #add_test_C_simple( ${TestName} )
ENDFOREACH( TestName  ${TestList} )



