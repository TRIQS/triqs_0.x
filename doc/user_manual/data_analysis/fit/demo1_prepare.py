import numpy as np
from pytriqs.Base.Plot.MatplotlibInterface import oplot
from pytriqs.Base.GF_Local import *
from pytriqs.Base.GF_Local.Descriptors import iOmega_n
g = GFBloc_ImFreq(Indices = [1], Beta = 300, NFreqMatsubara = 1000, Name = "g")

from pytriqs.Base.Archive import HDF_Archive
R = HDF_Archive('myfile.h5', 'w')

for n, Z0 in enumerate( np.arange (1,0, -0.1) ) :
    g <<= inverse( iOmega_n + 0.5  - iOmega_n * ( 1 - 1/Z0) ) # / (1 + 4*iOmega_n*iOmega_n) ) 
    g.Name = "Z = %s"%Z0 
    R[ str(n) ] = { 'Z0' : Z0, 'g' : g} 

