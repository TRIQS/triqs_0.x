from pytriqs.base.archive import HDFArchive
from pytriqs.base.gf_local import GFBloc_ImFreq

R = HDFArchive('myfile.h5', 'r')  # Opens the file myfile.h5 in readonly mode 
G = R['g1'] # Retrieve the object named g1 in the file as G

# ... ok now I can work with G


