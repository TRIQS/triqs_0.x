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

        GfImTime(Indices, Beta, Statistic, NTimeSlices,  Data, Tail, Name)

              * ``Indices``:  a list of indices names of the block
              * ``Beta``:  Inverse Temperature 
              * ``Statistic``:  GF_Statistic.Fermion [default] or GF_Statistic.Boson
              * ``NTimeTimes``  : Number of time points in the mesh
              * ``Data``:   A numpy array of dimensions (len(Indices),len(Indices),NTimeSlices) representing the value of the Green function on the mesh. 
              * ``Tail``:  the tail 
              * ``Name``:  a name of the GF
        If you already have the mesh, you can use a simpler version :

        GfImTime (Indices, Mesh, Data, Tail, Name)
           
              * ``Indices``:  a list of indices names of the block
              * ``Mesh``:  a MeshGf object, such that Mesh.TypeGF== GF_Type.Imaginary_Time 
              * ``Data``:   A numpy array of dimensions (len(Indices),len(Indices),NTimeSlices) representing the value of the Green function on the mesh. 
              * ``Tail``:  the tail 
              * ``Name``:  a name of the GF

        .. warning::
        The Green function take a **view** of the array Data, and a **reference** to the Tail.

        """
        mesh = d.pop('Mesh',None)
        if mesh is None : 
            if 'Beta' not in d : raise ValueError, "Beta not provided"
            Beta = float(d.pop('Beta'))
            Nmax = d.pop('NTimePoints',10000)
            stat = d.pop('Statistic','F') 
            sh = 1 if stat== 'F' else 0 
            mesh = MeshImTime(Beta,'F',Nmax)

        self.dtype = numpy.float64
        indicesL, indicesR = get_indices_in_dict(d)
        N1, N2 = len(indicesL),len(indicesR)
        data = d.pop('Data') if 'Data' in d else numpy.zeros((N1,N2,len(mesh)), self.dtype )
        tail= d.pop('Tail') if 'Tail' in d else TailGf( shape = (N1,N2), order_min=-1, size=10)
        symmetry = d.pop('Symmetry',None)
        name =  d.pop('Name','g')
        assert len(d) ==0, "Unknown parameters in GFBloc constructions %s"%d.keys() 
        
        GfImTime_cython.__init__(self, mesh, data, tail, symmetry,(indicesL,indicesR), name )

    #--------------   PLOT   ---------------------------------------
   
    def _plot_(self, OptionsDict):
        """ Plot protocol. OptionsDict can contain : 
             * :param RI: 'R', 'I', 'RI' [ default] 
             * :param x_window: (xmin,xmax) or None [default]
             * :param Name: a string [default ='']. If not '', it remplaces the name of the function just for this plot.
        """
        has_complex_value = False
        return impl_plot.plot_base( self, OptionsDict,  r'$\tau$', lambda name : r'%s$(\tau)$'%name, has_complex_value ,  list(self.mesh) )

