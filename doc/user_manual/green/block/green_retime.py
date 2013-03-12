import numpy
from pytriqs.gf.local import GfReFreq,SemiCircular
from pytriqs.plot.mpl_interface import oplot

gf = GfReFreq(indices = [0], omega_min = -30, omega_max = 30, n_freq_points = 1000, name = "my_block")
gf <<= SemiCircular(half_bandwidth = 2)

oplot(gf.imag, '-o')


