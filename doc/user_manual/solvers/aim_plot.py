from pytriqs.base.GF_Local import *
from pytriqs.base.Archive import *
from pytriqs.base.Plot.MatplotlibInterface import oplot

A = HDF_Archive("solution.h5")
oplot(A['G']['up'], '-o', x_window = (0,10))
