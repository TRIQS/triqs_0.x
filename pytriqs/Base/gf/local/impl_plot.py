from pytriqs.Base.Plot.protocol import clip_array

class _Plot_Wrapper_Partial_Reduce : 
    """ Internal Use"""
    def __init__(self, obj,  **opt) : 
        self.obj, self.opt = obj,opt
    def _plot_(self,Options) : 
        Options.update(self.opt)
        return self.obj._plot_(Options)

#----------------

def plot_base (self, OptionsDict, xlabel, ylabel, use_ris, X):
    """ Plot protocol. OptionsDict can contain : 
         * :param RIS: 'R', 'I', 'S', 'RI' [ default] 
         * :param x_window: (xmin,xmax) or None [default]
         * :param Name: a string [default ='']. If not '', it remplaces the name of the function just for this plot.
    """
    Name = OptionsDict.pop('Name', '' )  # consume it
    NamePrefix = OptionsDict.pop('NamePrefix', '' )  # consume it
    if Name and NamePrefix : raise ValueError, 'Name and NamePrefix can not be used at the same time'
    if NamePrefix : name_save, self.Name = self.Name, Name or NamePrefix

    rx = OptionsDict.pop('x_window',None ) # consume it
    sl = clip_array (X, *rx) if rx else slice(len(X)) # the slice due to clip option x_window

    def mdic( prefix, f) : 
       return [{'Type' : "XY", 
                'xlabel' : xlabel,
                'ylabel' : ylabel (self.Name),
                'xdata' : X[sl],
                'label' : Name if Name else prefix + B.Name ,
                'ydata' : f( B._data.array[0,0,sl] ) } for (i,j,B) in self ] 

    if use_ris : 
        ris = OptionsDict.pop('RI','RI') 
        if   ris == "R" : 
            res = mdic( 'Re ', lambda x : x.real)
        elif ris == "I" : 
            res = mdic( 'Im ', lambda x : x.imag)
        elif ris == "S" :
            res = mdic( '', lambda x : -1/numpy.pi *x.imag)
        elif ris == 'RI' :
             res = mdic( 'Re ', lambda x : x.real) + mdic( 'Im ', lambda x : x.imag)
        else : 
             raise ValueError, "RIS flags meaningless %s"%ris
    else: 
        res = mdic( '', lambda x : x)
        
    if NamePrefix: self.Name = name_save
    return res 

#------------------

def x_data_view (self, X, x_window = None, flatten_y = False) : 
    """
     :param x_window: the window of x variable (omega/omega_n/t/tau) for which data is requested
                      if None, take the full window
     :param flatten_y: If the Green function is of size (1,1) flatten the array as a 1d array
     :rtype: a tuple (X, data) where 
             * X is a 1d numpy of the x variable inside the window requested
             * data is a 3d numpy array of dim (:,:, len(X)), the corresponding slice of data
               If flatten_y is True and dim is (1,1,*), returns a 1d numpy
    """
    X = [x.imag for x in self._mesh] if self._mesh.TypeGF == GF_Type.Imaginary_Frequency  else [x for x in self._mesh]
    X, data = numpy.array(X), self._data.array
    if x_window :
      sl = clip_array (X, *x_window) if x_window else slice(len(X)) # the slice due to clip option x_window
      X, data = X[sl],  data[:,:,sl]
    if flatten_y and data.shape[:2]==(1,1) : data = data[0,0,:]
    return X, data
        
