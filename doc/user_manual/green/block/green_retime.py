from pytriqs.gf.local import GfReTime

g = GfReTime(indices = ['eg1', 'eg2'], window = (-5, 5), n_points = 1000, name = "egBlock")

#Example to be written with fourier transform!

#from pytriqs.plot.mpl_interface import oplot
#oplot(g['eg1','eg1'], '-o', RI = 'S')
#oplot(g['eg2','eg2'], '-x', RI = 'S')


