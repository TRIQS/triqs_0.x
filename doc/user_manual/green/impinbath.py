from pytriqs.base.gf_local import GFBloc_ReFreq, Omega, Wilson, inverse
import numpy
a = numpy.arange(-1.99,2.00,0.02) # Define the energy array
eps_d,V  = 0.3, 0.2

# Create the real-frequency Green's function and initialize it
g = GFBloc_ReFreq(Indices = ['s','d'], Beta = 50, MeshArray = a, Name = "s+d")
g['d','d'] = Omega - eps_d
g['d','s'] = V
g['s','d'] = V
g['s','s'] = inverse( Wilson(1.0) )
g.invert()

# Plot it with matplotlib. 'S' means: spectral function ( -1/pi Imag (g) )
from pytriqs.base.plot.mpl_interface import oplot
oplot( g['d','d'], '-o', RI = 'S', x_window  = (-1.8,1.8), Name = "Impurity" )
oplot( g['s','s'], '-x', RI = 'S', x_window  = (-1.8,1.8), Name = "Bath" )
