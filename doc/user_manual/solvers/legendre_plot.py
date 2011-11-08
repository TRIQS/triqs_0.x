from pytriqs.Base.GF_Local import *
from pytriqs.Base.Archive import *
from pytriqs.Base.Plot.MatplotlibInterface import oplot

A = HDF_Archive("solution.h5")
oplot(A['Gl']['up'], '-o', x_window=(15,45) )
