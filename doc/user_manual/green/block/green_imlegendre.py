from pytriqs.Base.GF_Local import *
from pytriqs.Base.Plot.MatplotlibInterface import oplot,plt

# A Green's function on the Matsubara axis set to a semicircular
gw = GFBloc_ImFreq(Indices = [1], Beta = 50)
gw <<= SemiCircular(HalfBandwidth = 1)

# Create a Legendre Green's function with 40 coefficients
# and initialize it from gw
gl = GFBloc_ImLegendre(Indices = [1], Beta = 50, NLegendreCoeffs = 40)
gl <<= MatsubaraToLegendre(gw)

# Plot the Legendre Green's function
oplot(gl, '-o')
