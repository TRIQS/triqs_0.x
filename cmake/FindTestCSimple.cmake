# runs python_with_DMFT script > output 
# and compares output with script.output
# Example: 
#   add_triqs_test_script(ExampleTest)
#   where ExampleTest.py is the script and ExampleTest.output is the expected output
#
macro(add_test_C_simple testname ) 
 enable_testing()
 add_test(${testname}${ARGN}
  ${CMAKE_COMMAND}
  -Dname=${testname}
  -Dcmd=${CMAKE_CURRENT_BINARY_DIR}/${testname}${ARGN}
  -Dreference=${CMAKE_CURRENT_SOURCE_DIR}/${testname}.output
  -P ${TRIQS_SOURCE_DIR}/cmake/run_test_simple.cmake
  )
endmacro(add_test_C_simple)


