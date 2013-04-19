from pytriqs.gf.local import *
from pytriqs.archive import *
from pytriqs.plot.mpl_interface import oplot

A = HDFArchive("solution.h5")
oplot(A['Gl']['up'], '-o', x_window=(15,45) )
