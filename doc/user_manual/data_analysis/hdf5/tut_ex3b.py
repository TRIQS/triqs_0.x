from pytriqs.Base.GF_Local import GFBloc_ReFreq
from pytriqs.Base.Archive import HDF_Archive
from math import pi

R = HDF_Archive('myfile.h5', 'r') 
 
from pytriqs.Base.Plot.MatplotlibPlotter import oplot, plt
plt.xrange(-1,1) 
plt.yrange(0,7) 

for name, g in R.items() :  # iterate on the elements of R, like a dict ...
    oplot( (- 1/pi * g).imag, "-o", Name = name)

p.savefig("./tut_ex3b.png") 

