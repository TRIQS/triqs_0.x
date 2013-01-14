###############  Density  #########################

cdef extern from "triqs/gf/local/density.hpp":
    matrix_view density(gf_imfreq &)

