from gf import GfImFreq_cython, MeshImFreq, TailGf
from gf_generic import GfGeneric
import numpy
from tools import get_indices_in_dict

class GfImFreq ( GfImFreq_cython, GfGeneric ) :
    def __init__(self, **d): 
        """
        The constructor have two variants : you can either provide the mesh in
        Matsubara frequencies yourself, or give the parameters to build it.
        All parameters must be given with keyword arguments.

        GfImFreq(Indices, Beta, Statistic, NFreqMatsubara,  Data, Tail, Name)

              * ``Indices``:  a list of indices names of the block
              * ``Beta``:  Inverse Temperature 
              * ``Statistic``:  GF_Statistic.Fermion [default] or GF_Statistic.Boson
              * ``NFreqMatsubara``:  Number of Matsubara frequencies
              * ``Data``:   A numpy array of dimensions (len(Indices),len(Indices),NFreqMatsubara) representing the value of the Green function on the mesh. 
              * ``Tail``:  the tail 
              * ``Name``:  a name of the GF

        If you already have the mesh, you can use a simpler version :

        GfImFreq(Indices, Mesh, Data, Tail, Name)
            
              * ``Indices``:  a list of indices names of the block
              * ``Mesh``:  a MeshGF object, such that Mesh.TypeGF== GF_Type.Imaginary_Frequency 
              * ``Data``:   A numpy array of dimensions (len(Indices),len(Indices),NFreqMatsubara) representing the value of the Green function on the mesh. 
              * ``Tail``:  the tail 
              * ``Name``:  a name of the GF

        .. warning::
        The Green function take a **view** of the array Data, and a **reference** to the Tail.
        """
        mesh = d.pop('Mesh',None)
        if mesh is None : 
            if 'Beta' not in d : raise ValueError, "Beta not provided"
            Beta = float(d.pop('Beta'))
            Nmax = d.pop('NFreqMatsubara',1025)
            stat = d.pop('Statistic','F') 
            sh = 1 if stat== 'F' else 0 
            mesh = MeshImFreq(Beta,'F',Nmax)

        self.dtype = numpy.complex_
        indicesL, indicesR = get_indices_in_dict(d)
        N1, N2 = len(indicesL),len(indicesR)
        data = d.pop('Data') if 'Data' in d else numpy.zeros((N1,N2,len(mesh)), self.dtype )
        tail= d.pop('Tail') if 'Tail' in d else TailGf( shape = (N1,N2), order_min=-1, size=10)
        symmetry = d.pop('Symmetry',None)
        name =  d.pop('Name','g')
        assert len(d) ==0, "Unknown parameters in GFBloc constructions %s"%d.keys() 
        
        GfImFreq_cython.__init__(self, mesh, data, tail, symmetry,(indicesL,indicesR), name )
    
    #--------------   PLOT   ---------------------------------------
   
    def _plot_(self, OptionsDict):
        """ Plot protocol. OptionsDict can contain : 
             * :param RIS: 'R', 'I', 'S', 'RI' [ default] 
             * :param x_window: (xmin,xmax) or None [default]
             * :param Name: a string [default ='']. If not '', it remplaces the name of the function just for this plot.
        """
        return impl_plot.plot_base( self, OptionsDict,  r'$\omega_n$', 
                lambda name : r'%s$(i\omega_n)$'%name, True, [x.imag for x in self.mesh] )
    #--------------   OTHER OPERATIONS -----------------------------------------------------

    def replace_by_tail(self,start) : 
        d = self.data
        t = self.tail
        for n, om in enumerate(self.mesh) : 
            if n >= start : d[:,:,n] = t(om).array

