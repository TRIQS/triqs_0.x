from pytriqs.lattice.tight_binding import *
from pytriqs.dos import HilbertTransform
from pytriqs.gf.local import GfImFreq

# Define a DOS (here on a square lattice)
BL = BravaisLattice(Units = [(1,0,0) , (0,1,0) ], orbital_positions= {"" :  (0,0,0)} ) 
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

TB = TightBinding (BL, hop)
d = dos(TB, n_kpts= 500, n_eps = 101, name = 'dos')[0]

#define a Hilbert transform
H = HilbertTransform(d)

#fill a Green function
G = GfImFreq(indices = ['up','down'], beta = 20)
Sigma0 = GfImFreq(indices = ['up','down'], beta = 20); Sigma0.zero()
G <<= H(Sigma = Sigma0,mu=0.)

