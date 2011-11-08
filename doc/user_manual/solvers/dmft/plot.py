from pytriqs.Base.GF_Local import *
from pytriqs.Base.Archive import *
from pytriqs.Base.Plot.MatplotlibInterface import *

A = HDF_Archive("SingleSiteBethe.h5",'r')

for i in range(5):
  oplot(A['G-%s'%i]['up'].imag,'-o',Name='Iteration = %s'%i, x_window = (0,2))

