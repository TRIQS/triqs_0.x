SET(PYTHON_SOURCES 
 ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/models.py 
)

# Install python sources
install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/applications/impurity_solvers/ctqmc_hyb)

# Build cython module 
cython_module(ctqmc_hyb ctqmc_solver applications/impurity_solvers/ctqmc_hyb)
target_link_libraries(ctqmc_solver triqs_solver_ctqmc)


