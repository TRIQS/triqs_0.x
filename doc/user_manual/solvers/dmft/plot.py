from pytriqs.gf.local import *
from pytriqs.archive import *
from pytriqs.plot.mpl_interface import *

A = HDFArchive("single_site_bethe.h5",'r')

for i in range(5):
  oplot(A['G-%s'%i]['up'].imag,'-o', name='Iteration = %s'%i, x_window = (0,2))

