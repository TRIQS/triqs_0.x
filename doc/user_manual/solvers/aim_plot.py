from pytriqs.base.gf.local import *
from pytriqs.base.archive import *
from pytriqs.base.plot.mpl_interface import oplot

A = HDFArchive("solution.h5")
oplot(A['G']['up'], '-o', x_window = (0,10))
