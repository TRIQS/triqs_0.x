import numpy as np
from pytriqs.base.GF_Local import GFBloc_ReFreq, SemiCircular

g = GFBloc_ReFreq(Indices = ['eg1', 'eg2'], Beta = 50, MeshArray = np.arange(-5,5,0.01) , Name = "egBlock")

g['eg1','eg1'] = SemiCircular(HalfBandwidth = 1)
g['eg2','eg2'] = SemiCircular(HalfBandwidth = 2)

from pytriqs.base.Plot.MatplotlibInterface import oplot
oplot(g['eg1','eg1'], '-o', RI = 'S')  # S : spectral function 
oplot(g['eg2','eg2'], '-x', RI = 'S')   


