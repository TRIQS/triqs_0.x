from gf import GfImFreq_cython, MeshImFreq, TailGf
from gf_generic import GfGeneric
import numpy
from tools import get_indices_in_dict
import impl_plot

class GfImFreq ( GfImFreq_cython, GfGeneric ) :
    def __init__(self, **d): 
        """
        The constructor have two variants : you can either provide the mesh in
        Matsubara frequencies yourself, or give the parameters to build it.
        All parameters must be given with keyword arguments.

        GfImFreq(indices, beta, statistic, n_matsubara, data, tail, name)

              * ``indices``:  a list of indices names of the block
              * ``beta``:  Inverse Temperature 
              * ``statistic``:  'F' or 'B'
              * ``n_matsubara``:  Number of Matsubara frequencies
              * ``data``:   A numpy array of dimensions (len(indices),len(indices),n_matsubara) representing the value of the Green function on the mesh. 
              * ``tail``:  the tail 
              * ``name``:  a name of the GF

        If you already have the mesh, you can use a simpler version :

        GfImFreq(indices, mesh, data, tail, name)
            
              * ``indices``:  a list of indices names of the block
              * ``mesh``:  a MeshGf object, such that mesh.TypeGF== GF_Type.Imaginary_Frequency 
              * ``data``:   A numpy array of dimensions (len(indices),len(indices),n_matsubara) representing the value of the Green function on the mesh. 
              * ``tail``:  the tail 
              * ``name``:  a name of the GF

        .. warning::
        The Green function take a **view** of the array data, and a **reference** to the tail.
        """
        mesh = d.pop('mesh',None)
        if mesh is None : 
            if 'beta' not in d : raise ValueError, "beta not provided"
            beta = float(d.pop('beta'))
            n_max = d.pop('n_matsubara',1025)
            stat = d.pop('statistic','F') 
            sh = 1 if stat== 'F' else 0 
            mesh = MeshImFreq(beta,'F',n_max)

        self.dtype = numpy.complex_
        indicesL, indicesR = get_indices_in_dict(d)
        N1, N2 = len(indicesL),len(indicesR)
        data = d.pop('data') if 'data' in d else numpy.zeros((N1,N2,len(mesh)), self.dtype )
        tail= d.pop('tail') if 'tail' in d else TailGf(shape = (N1,N2), size=10, order_min=-1)
        symmetry = d.pop('symmetry', None)
        name =  d.pop('name','g')
        assert len(d) ==0, "Unknown parameters in GFBloc constructions %s"%d.keys() 
        
        GfImFreq_cython.__init__(self, mesh, data, tail, symmetry,(indicesL,indicesR), name )
    
    #--------------   PLOT   ---------------------------------------
   
    def _plot_(self, opt_dict):
        """ Plot protocol. opt_dict can contain : 
             * :param RIS: 'R', 'I', 'S', 'RI' [ default] 
             * :param x_window: (xmin,xmax) or None [default]
             * :param name: a string [default ='']. If not '', it remplaces the name of the function just for this plot.
        """
        return impl_plot.plot_base( self, opt_dict,  r'$\omega_n$', 
                lambda name : r'%s$(i\omega_n)$'%name, True, [x.imag for x in self.mesh] )
    #--------------   OTHER OPERATIONS -----------------------------------------------------

    def replace_by_tail(self,start) : 
        d = self.data
        t = self.tail
        for n, om in enumerate(self.mesh) : 
            if n >= start : d[:,:,n] = t(om).array

