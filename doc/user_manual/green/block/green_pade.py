import numpy
from math import pi
from cmath import sqrt, log
from pytriqs.base.gf.local import *
from pytriqs.base.gf.local.descriptors import Function

beta = 100  # Inverse temperature
L = 101     # Number of Matsubara frequencies used in the Pade approximation
eta = 0.01  # Imaginary frequency shift

## Test Green's functions ##

# Two Lorentzians
def GLorentz(z):
    return 0.7/(z-2.6+0.3*1j) + 0.3/(z+3.4+0.1*1j)

# Semicircle
def GSC(z):
    return 2.0*(z + sqrt(1-z**2)*(log(1-z) - log(-1+z))/pi)

# A superposition of GLorentz(z) and GSC(z) with equal weights
def G(z):
    return 0.5*GLorentz(z) + 0.5*GSC(z)

# Matsubara GF
gm = GfImFreq(indices = [0], beta = beta, name = "gm")
gm <<= Function(G)
gm.tail.zero()
gm.tail[1] = numpy.array([[1.0]])

# Real frequency BlockGf(reference)
gr = GfReFreq(indices = [0], omega_min = -6, omega_max = 6, n_freq_points = 1000, name = "gr")
gr <<= Function(G)
gr.tail.zero()
gr.tail[1] = numpy.array([[1.0]])

# Analytic continuation of gm
g_pade = GfReFreq(indices = [0], omega_min = -6, omega_max = 6, n_freq_points = 1000, name = "g_pade")
g_pade.set_from_pade(gm, N_Matsubara_Frequencies = L, Freq_Offset = eta)

# Comparison plot
from pytriqs.base.plot.mpl_interface import oplot
oplot(gr[0,0], '-o', RI = 'S', name = "Original DOS")
oplot(g_pade[0,0], '-x', RI = 'S', name = "Pade-reconstructed DOS")
