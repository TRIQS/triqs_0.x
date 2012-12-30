###############  Fourier  #########################

cdef extern from "triqs/gf/local/fourier_matsubara.hpp" : 
    gf_imfreq lazy_fourier          (gf_imtime & )
    gf_imtime lazy_inverse_fourier  (gf_imfreq & )


