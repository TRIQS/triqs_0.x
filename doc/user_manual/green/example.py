# Import the Green's functions 
from pytriqs.base.gf_local import GFBloc_ImFreq, iOmega_n, inverse 

# Create the Matsubara-frequency Green's function and initialize it
g = GFBloc_ImFreq(Indices = [1], Beta = 50, NFreqMatsubara = 1000, Name = "imp")
g <<= inverse( iOmega_n + 0.5 )

from pytriqs.base.plot.mpl_interface import oplot
oplot(g, '-o',  x_window  = (0,10))

