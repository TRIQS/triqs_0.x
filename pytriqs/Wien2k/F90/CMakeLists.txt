 find_package(BuildF2pyModule)
 BuildF2pyModule( triqs_Wien2k vertex vertex.pyf vertex.f90)
 install (FILES ${CMAKE_CURRENT_BINARY_DIR}/vertex.so DESTINATION ${TRIQS_PYTHON_LIB_DEST}/Wien2k)
 execute_process(COMMAND ln -fs ${CMAKE_CURRENT_BINARY_DIR}/vertex.so  ${CMAKE_CURRENT_BINARY_DIR}/.. )
