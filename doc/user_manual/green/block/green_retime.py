import numpy
from pytriqs.base.GF_Local import GFBloc_ReFreq,SemiCircular
from pytriqs.base.Plot.MatplotlibInterface import oplot

gf = GFBloc_ReFreq(Indices = [1], Beta = 50, MeshArray = numpy.arange(-30,30,0.1) , Name = "my_block")
gf[1,1] =  SemiCircular(HalfBandwidth = 2)

oplot(gf.InverseFourier().imag, '-o')   


