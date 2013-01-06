import numpy
from pytriqs.base.gf_local import GfReFreq,SemiCircular
from pytriqs.base.plot.mpl_interface import oplot

gf = GfReFreq(indices = [1], beta = 50, mesh_array = numpy.arange(-30,30,0.1) , name = "my_block")
gf[1,1] =  SemiCircular(half_bandwidth = 2)

oplot(gf.InverseFourier().imag, '-o')   


