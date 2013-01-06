import numpy as np
from pytriqs.base.gf_local import GfReFreq, SemiCircular

g = GfReFreq(indices = ['eg1', 'eg2'], beta = 50, mesh_array = np.arange(-5,5,0.01) , name = "egBlock")

g['eg1','eg1'] = SemiCircular(half_bandwidth = 1)
g['eg2','eg2'] = SemiCircular(half_bandwidth = 2)

from pytriqs.base.plot.mpl_interface import oplot
oplot(g['eg1','eg1'], '-o', RI = 'S')  # S : spectral function 
oplot(g['eg2','eg2'], '-x', RI = 'S')   


