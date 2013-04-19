# runs python_with_DMFT script > output 
# and compares output with script.output
# Example: 
#   add_triqs_test_script(ExampleTest)
#   where ExampleTest.py is the script and ExampleTest.output is the expected output
#
macro(add_test_C_simple testname ) 
 enable_testing()

 if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${testname}.output)

  file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/${testname}.output DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  add_test(${testname}
   ${CMAKE_COMMAND}
   -Dname=${testname}${ARGN}
   -Dcmd=${CMAKE_CURRENT_BINARY_DIR}/${testname}${ARGN}
   -Dreference=${CMAKE_CURRENT_SOURCE_DIR}/${testname}.output
   -P ${TRIQS_SOURCE_DIR}/cmake/run_test_simple.cmake
   )

 else (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${testname}.output)
  add_test(${testname}${ARGN} ${testname}${ARGN} )
 endif (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${testname}.output)

endmacro(add_test_C_simple)


