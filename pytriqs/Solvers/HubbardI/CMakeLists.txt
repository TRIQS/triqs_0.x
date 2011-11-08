set(NickName "HubbardI")

# install the python sources 
SET(PYTHON_SOURCES 
 ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/Solver_HubbardI.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/Solver_HubbardI_base.py 
 )

# 
python_register_dynamic_module(${NickName})

# build the fortran stuff...
if (Build_${NickName})
  install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/Solvers/HubbardI)
  add_subdirectory(F90)
  add_subdirectory(test)
endif (Build_${NickName})
