from pytriqs.Base.Archive import HDF_Archive
from pytriqs.Base.GF_Local import GFBloc_ImFreq

R = HDF_Archive('myfile.h5', 'r')  # Opens the file myfile.h5 in readonly mode 
G = R['g1'] # Retrieve the object named g1 in the file as G

# ... ok now I can work with G


