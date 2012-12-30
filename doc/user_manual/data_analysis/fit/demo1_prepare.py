import numpy as np
from pytriqs.base.plot.mpl_interface import oplot
from pytriqs.base.gf_local import *
from pytriqs.base.gf_local.descriptors import iOmega_n
g = GFBloc_ImFreq(Indices = [1], Beta = 300, NFreqMatsubara = 1000, Name = "g")

from pytriqs.base.archive import HDF_Archive
R = HDF_Archive('myfile.h5', 'w')

for n, Z0 in enumerate( np.arange (1,0, -0.1) ) :
    g <<= inverse( iOmega_n + 0.5  - iOmega_n * ( 1 - 1/Z0) ) # / (1 + 4*iOmega_n*iOmega_n) ) 
    g.Name = "Z = %s"%Z0 
    R[ str(n) ] = { 'Z0' : Z0, 'g' : g} 

