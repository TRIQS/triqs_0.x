from gf import GfImTime_cython, MeshImTime, TailGf
from gf_generic import GfGeneric
import numpy
from tools import get_indices_in_dict
import impl_plot

class GfImTime ( GfImTime_cython, GfGeneric ) :
    def __init__(self, **d):
        """
        The constructor have two variants : you can either provide the mesh in
        Matsubara frequencies yourself, or give the parameters to build it.
        All parameters must be given with keyword arguments.

        GfImTime(indices, beta, statistic, n_time_points, data, tail, name)

              * ``indices``:  a list of indices names of the block
              * ``beta``:  Inverse Temperature 
              * ``statistic``:  'F' or 'B'
              * ``n_time_points``  : Number of time points in the mesh
              * ``data``:   A numpy array of dimensions (len(indices),len(indices),n_time_points) representing the value of the Green function on the mesh. 
              * ``tail``:  the tail 
              * ``name``:  a name of the GF

        GfImTime (indices, mesh, data, tail, name)
           
              * ``indices``:  a list of indices names of the block
              * ``mesh``:  a MeshGf object, such that mesh.TypeGF== GF_Type.Imaginary_Time 
              * ``data``:   A numpy array of dimensions (len(indices),len(indices),n_time_points) representing the value of the Green function on the mesh. 
              * ``tail``:  the tail 
              * ``name``:  a name of the GF

        .. warning::

          The Green function take a **view** of the array data, and a **reference** to the tail.

        """
        mesh = d.pop('mesh',None)
        if mesh is None : 
            if 'beta' not in d : raise ValueError, "beta not provided"
            beta = float(d.pop('beta'))
            stat = d.pop('statistic','F') 
            n_max = d.pop('n_time_points',10000)
            kind = d.pop('kind','H') 
            mesh = MeshImTime(beta,stat,n_max,kind)

        self.dtype = numpy.float64
        indicesL, indicesR = get_indices_in_dict(d)
        N1, N2 = len(indicesL),len(indicesR)
        data = d.pop('data') if 'data' in d else numpy.zeros((N1,N2,len(mesh)), self.dtype )
        tail= d.pop('tail') if 'tail' in d else TailGf(shape = (N1,N2), size=10,  order_min=-1)
        symmetry = d.pop('symmetry',None)
        name =  d.pop('name','g')
        assert len(d) ==0, "Unknown parameters in GFBloc constructions %s"%d.keys() 
        
        GfImTime_cython.__init__(self, mesh, data, tail, symmetry, (indicesL,indicesR), name)

    #--------------   PLOT   ---------------------------------------
   
    def _plot_(self, OptionsDict):
        """ Plot protocol. OptionsDict can contain : 
             * :param RI: 'R', 'I', 'RI' [ default] 
             * :param x_window: (xmin,xmax) or None [default]
             * :param name: a string [default ='']. If not '', it remplaces the name of the function just for this plot.
        """
        has_complex_value = False
        return impl_plot.plot_base( self, OptionsDict,  r'$\tau$', lambda name : r'%s$(\tau)$'%name, has_complex_value ,  list(self.mesh) )

