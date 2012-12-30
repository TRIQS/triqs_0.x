from pytriqs.base.gf_local import *
from pytriqs.base.gf_local.descriptors import Omega
g = GFBloc_ImFreq(Indices = [1], Beta = 50, NFreqMatsubara = 1000, Name = "g")
g <<= inverse( Omega + 0.5 )

# open 2 panels top (t) and bottom (b) 
from pytriqs.base.plot.mpl_interface import subplots
f, (t,b) = subplots( 2,1)

#plot ...
t.oplot(g.real, '-o', x_window = (0,10) )
b.oplot(g.imag, '-x', x_window = (0,12) )   
b.oplot( lambda om : -om*0.8/(om*om + 4), Name = "Bad Fit !")
b.text(5,-0.5, r'$g(i\omega_n) = \frac{1}{i \omega_n + 1/2} $', size = 20, color='r')

