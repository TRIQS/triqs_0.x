from pytriqs.base.gf.local import *
from pytriqs.base.plot.mpl_interface import oplot,plt

# A Green's function on the Matsubara axis set to a semicircular
gw = GfImFreq(indices = [1], beta = 50)
gw <<= SemiCircular(half_bandwidth = 1)

# Create an imaginary-time Green's function
gt = GfImTime(indices = [1], beta = 50)
gt <<= InverseFourier(gw)

# Plot the Legendre Green's function
oplot(gt, '-')
