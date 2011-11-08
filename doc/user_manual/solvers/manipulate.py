from pytriqs.Solvers.Operators import *

H = C('up',1) * Cdag('up',2) + C('up',2) * Cdag('up',1)
print H
print H - H.dagger()

print AntiCommutator(C('up'),Cdag('up'))
print AntiCommutator(C('up'),0.5*Cdag('down'))
