from pytriqs.base.lattice.tight_binding import *
from pytriqs.base.dos import Hilbert_Transform
from pytriqs.base.gf_local import GfImFreq

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
d = dos (TB, nkpts= 500, neps = 101, name = 'dos')[0]

#define a Hilbert transform
H = Hilbert_Transform(d)

#fill a Green function
G = GfImFreq(indices = ['up','down'], beta = 20)
Sigma0 = GfImFreq(indices = ['up','down'], beta = 20); Sigma0.zero()
G <<= H(Sigma = Sigma0,mu=0.)

