find_package(BuildF2pyModule)
BuildF2pyModule( triqs_HubbardI hubbard_I hubbard_I.pyf hubbard_I.f90)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/hubbard_I.so DESTINATION ${TRIQS_PYTHON_LIB_DEST}/applications/impurity_solvers/hubbard_I)

