set(S triqs_solver_ctqmc)

add_library(${S} configuration.cpp dynamic_trace.cpp hloc.cpp ctqmc.cpp measure_legendre.cpp measure_correlator.cpp)
target_link_libraries (${S} triqs)

install(TARGETS ${S} DESTINATION lib)

