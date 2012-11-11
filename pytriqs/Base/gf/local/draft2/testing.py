import numpy as np
import _pytriqs_GF2
import TailGF

#T = TailGF.TailGF(IndicesL=[1],IndicesR=[1], size=3, OrderMin=-1)
T = TailGF.TailGF(IndicesL=[1],IndicesR=[1], array=np.array([[[2,3,4,5,6]]],np.complex,order='F'), OrderMin=-1)
print T

M = _pytriqs_GF2.MeshMatsubaraFrequency(10.0,_pytriqs_GF2.GF_Statistic.Fermion,1025)

G = _pytriqs_GF2.GF(M,np.zeros((1,1,10),np.complex,order='F'),T)

