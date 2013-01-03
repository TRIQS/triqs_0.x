from pytriqs.base.gf_local import GFBloc_ReFreq
from pytriqs.base.archive import HDFArchive
from math import pi

R = HDFArchive('myfile.h5', 'r') 
 
from pytriqs.base.plot.mpl_interface import oplot, plt
plt.xrange(-1,1) 
plt.yrange(0,7) 

for name, g in R.items() :  # iterate on the elements of R, like a dict ...
    oplot( (- 1/pi * g).imag, "-o", Name = name)

p.savefig("./tut_ex3b.png") 

