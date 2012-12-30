 find_package(BuildF2pyModule)
 BuildF2pyModule( triqs_DFT vertex vertex.pyf vertex.f90)
 install (FILES ${CMAKE_CURRENT_BINARY_DIR}/vertex.so DESTINATION ${TRIQS_PYTHON_LIB_DEST}/dft)
 execute_process(COMMAND ln -fs ${CMAKE_CURRENT_BINARY_DIR}/vertex.so  ${CMAKE_CURRENT_BINARY_DIR}/.. )
