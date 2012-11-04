from arrays cimport *   
from shared_ptr cimport *
from dcomplex cimport * 
from gf_basic_tools cimport *
from tail cimport *
from gf_matsubara_freq cimport *
from gf_matsubara_time cimport *

cdef extern from "triqs/gf/local/fourier_matsubara.hpp" : 
    gf_view_freq lazy_fourier (gf_view_time & )
    gf_view_time lazy_inverse_fourier (gf_view_freq & )
    #gf_view[matsubara_freq_desc]  lazy_fourier          (gf_view_time & )
    #gf_view_t
