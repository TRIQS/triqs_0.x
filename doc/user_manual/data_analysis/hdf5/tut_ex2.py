from pytriqs.base.archive import HDF_Archive
from pytriqs.base.gf_local import GFBloc_ImFreq

# Define a Green function 
G = GFBloc_ImFreq ( Indices = [1], Beta = 10, NFreqMatsubara = 1000) 
      
# Opens the file myfile.h5, in read/write mode
R = HDF_Archive('myfile.h5', 'w')
# Store the object G under the name 'g1' and mu
R['g1'] = G
R['mu'] = 1.29
del R # closing the files (optional : file is closed when the R reference is deleted)




