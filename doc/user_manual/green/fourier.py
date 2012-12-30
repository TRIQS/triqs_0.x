from pytriqs.base.GF_Local import *
from pytriqs.base.Plot.MatplotlibInterface import oplot,plt

# A Green's function on the Matsubara axis set to a semicircular
gw = GFBloc_ImFreq(Indices = [1], Beta = 50)
gw <<= SemiCircular(HalfBandwidth = 1)

# Create an imaginary-time Green's function and plot it
gt = GFBloc_ImTime(Indices = [1], Beta = 50)
gt <<= InverseFourier(gw)
oplot(gt, '-')
