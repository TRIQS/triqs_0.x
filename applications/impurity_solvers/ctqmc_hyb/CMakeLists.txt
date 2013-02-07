set(S triqs_solver_ctqmc)
add_library(${S} Configuration.cpp DynamicTrace.cpp Hloc.cpp MC.cpp Measures_Legendre_allseries.cpp Measures_Legendre.cpp Measures_OpCorr.cpp)
target_link_libraries (${S} triqs)

#install (FILES ${CMAKE_CURRENT_BINARY_DIR}/ctqmc_hyb.so DESTINATION ${TRIQS_PYTHON_LIB_DEST}/solvers/ctqmc_hyb)
