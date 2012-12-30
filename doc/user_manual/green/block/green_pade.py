import numpy
from math import pi
from cmath import sqrt, log
from pytriqs.base.GF_Local import *
from pytriqs.base.GF_Local.Descriptors import Function

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
gm = GFBloc_ImFreq(Indices = [0], Beta = beta, Name = "gm")
gm <<= Function(G)
gm._tail.zero()
gm._tail[1] = numpy.array([[1.0]])

# Real frequency GF (reference)
gr = GFBloc_ReFreq(Indices = [0], Beta = beta, MeshArray = numpy.arange(-6,6,0.01), Name = "gr")
gr <<= Function(G)
gr._tail.zero()
gr._tail[1] = numpy.array([[1.0]])

# Analytic continuation of gm
g_pade = GFBloc_ReFreq(Indices = [0], Beta = beta, MeshArray = numpy.arange(-6,6,0.01), Name = "g_pade")
g_pade.setFromPadeOf(gm, N_Matsubara_Frequencies = L, Freq_Offset = eta)

# Comparison plot
from pytriqs.base.plot.MatplotlibInterface import oplot
oplot(gr[0,0], '-o', RI = 'S', Name = "Original DOS")
oplot(g_pade[0,0], '-x', RI = 'S', Name = "Pade-reconstructed DOS")
