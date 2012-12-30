from pytriqs.base.GF_Local import *
from pytriqs.base.plot.MatplotlibInterface import oplot,plt

# A Green's function on the Matsubara axis set to a semicircular
gw = GFBloc_ImFreq(Indices = [1], Beta = 50)
gw <<= SemiCircular(HalfBandwidth = 1)

# Create a Legendre Green's function with 40 coefficients
# initialize it from gw and plot it
gl = GFBloc_ImLegendre(Indices = [1], Beta = 50, NLegendreCoeffs = 40)
gl <<= MatsubaraToLegendre(gw)
oplot(gl, '-o')
