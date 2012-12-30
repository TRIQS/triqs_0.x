from pytriqs.base.lattice.TightBinding import *
from pytriqs.base.dos.Hilbert_Transform import *
from pytriqs.base.GF_Local import GFBloc_ImFreq

# Define a DOS (here on a square lattice)
BL = bravais_lattice(Units = [(1,0,0) , (0,1,0) ], Orbital_Positions= {"" :  (0,0,0)} ) 
t   = -1.00                # First neighbour Hopping
tp  =  0.0*t               # Second neighbour Hopping
hop= {  (1,0)  :  [[ t]],       
        (-1,0) :  [[ t]],     
        (0,1)  :  [[ t]],
        (0,-1) :  [[ t]],
        (1,1)  :  [[ tp]],
        (-1,-1):  [[ tp]],
        (1,-1) :  [[ tp]],
        (-1,1) :  [[ tp]]}
TB = tight_binding ( BL, hop)
d = dos (TB, nkpts= 500, neps = 101, Name = 'dos')[0]

#define a Hilbert transform
H = Hilbert_Transform(d)

#fill a Green function
G = GFBloc_ImFreq(Indices = ['up','down'], Beta = 20)
Sigma0 = GFBloc_ImFreq(Indices = ['up','down'], Beta = 20); Sigma0.zero()
G <<= H(Sigma = Sigma0,mu=0.)

