from pytriqs.Base.GF_Local import *
from pytriqs.Base.Archive import *
from pytriqs.Base.Plot.MatplotlibInterface import oplot

A = HDF_Archive("solution.h5")
oplot(A['G']['up'], '-o', x_window = (0,10))
