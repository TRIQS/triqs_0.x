set(NickName "HubbardI")

# install the python sources 
SET(PYTHON_SOURCES 
 ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/solver.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/solver_base.py 
 )

# 
python_register_dynamic_module(${NickName})

# build the fortran stuff...
if (Build_${NickName})
  install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/solvers/hubbard_I)
  add_subdirectory(test)
endif (Build_${NickName})
