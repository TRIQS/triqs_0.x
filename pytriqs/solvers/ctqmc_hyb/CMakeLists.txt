SET(PYTHON_SOURCES 
 ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/solver.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/models.py 
 )

option( Build_CTHyb "Should I build CTHyb ?" ON)
if (Build_CTHyb)
 install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/solvers/ctqmc_hyb)
 add_subdirectory(test)
 execute_process(COMMAND ln -fs ${CMAKE_BINARY_DIR}/applications/impurity_solvers/ctqmc_hyb/ctqmc_hyb.so ${CMAKE_CURRENT_BINARY_DIR})
endif (Build_CTHyb)

