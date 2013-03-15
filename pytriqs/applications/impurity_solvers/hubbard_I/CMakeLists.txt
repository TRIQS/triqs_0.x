set(NickName "HubbardI")

# install the python sources 
SET(PYTHON_SOURCES 
 ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/hubbard_solver.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/old_solver.py 
 )

# 
python_register_dynamic_module(${NickName})

# build the fortran stuff...
install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/applications/impurity_solvers/hubbard_I)
execute_process(COMMAND ln -fs ${CMAKE_BINARY_DIR}/applications/impurity_solvers/hubbard_I/hubbard_I.so ${CMAKE_CURRENT_BINARY_DIR})
