from pytriqs.base.gf_local import GFBloc_ImFreq, SemiCircular

g = GFBloc_ImFreq(Indices = ['eg1','eg2'], Beta = 50, NFreqMatsubara = 1000, Name = "egBlock") 

g['eg1','eg1'] = SemiCircular(HalfBandwidth = 1)
g['eg2','eg2'] = SemiCircular(HalfBandwidth = 2)

from pytriqs.base.plot.mpl_interface import oplot,plt
oplot(g, '-o', x_window = (0,10))
plt.ylim(-2,1)

