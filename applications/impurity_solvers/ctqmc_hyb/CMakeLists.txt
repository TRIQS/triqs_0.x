set( S ctqmc_hyb)
add_library(${S} Configuration.cpp DynamicTrace.cpp Hloc.cpp MC.cpp Measures_Legendre_allseries.cpp Measures_Legendre.cpp Measures_OpCorr.cpp QMC_wrap.cpp )
set_target_properties(${S} PROPERTIES PREFIX "")
target_link_libraries ( ${S}  triqs )

install (FILES ${CMAKE_CURRENT_BINARY_DIR}/ctqmc_hyb.so DESTINATION ${TRIQS_PYTHON_LIB_DEST}/solvers/ctqmc_hyb)

