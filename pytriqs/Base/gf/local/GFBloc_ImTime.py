################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Ferrero, O. Parcollet
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
__all__ = ['GFBloc_ImTime']
import numpy
from math import pi
from pytriqs_GF3 import GFBloc_ImTime_cython, MeshMatsubaraTime  
from GFBloc_general import _GFBloc_general 

class GFBloc_ImTime (GFBloc_ImTime_cython,  _GFBloc_general):
    """ 
    A matrix-valued block Green's function in Matsubara time.
    """
    def __init__(self, **d):
        """
  
     The constructor have two variants : you can either provide the mesh in
     Matsubara frequencies yourself, or give the parameters to build it.
     All parameters must be given with keyword arguments.

     GFBloc_ImTime(Indices, Beta, Statistic, NTimeSlices,  Data, Tail, Name,Note)
           * ``Indices``:  a list of indices names of the block
           * ``Beta``:  Inverse Temperature 
           * ``Statistic``:  GF_Statistic.Fermion [default] or GF_Statistic.Boson
           * ``NTimeSlices``  : Number of time slices in the mesh
           * ``Data``:   A numpy array of dimensions (len(Indices),len(Indices),NTimeSlices) representing the value of the Green function on the mesh. 
           * ``Tail``:  the tail 
           * ``Name``:  a name of the GF
           * ``Note``:  any string you like...

     If you already have the mesh, you can use a simpler version :

     GFBloc_ImTime (Indices, Mesh, Data, Tail, Name,Note)
        
           * ``Indices``:  a list of indices names of the block
           * ``Mesh``:  a MeshGF object, such that Mesh.TypeGF== GF_Type.Imaginary_Time 
           * ``Data``:   A numpy array of dimensions (len(Indices),len(Indices),NTimeSlices) representing the value of the Green function on the mesh. 
           * ``Tail``:  the tail 
           * ``Name``:  a name of the GF
           * ``Note``:  any string you like...

 .. warning::
    The Green function take a **view** of the array Data, and a **reference** to the Tail.

         """
        # construct the mesh if needed
        if 'Mesh' not in d : 
            if 'Beta' not in d : raise ValueError, "Beta not provided"
            Beta = float(d.pop('Beta'))
            Nmax = d.pop('NTimeSlices', 10000)
            stat = d.pop('Statistic','F')
            d['Mesh'] = MeshMatsubaraTime(Beta,'F',Nmax)
            # TO RECHECK 
            #d['Mesh'] = MeshGF( GF_Type.Imaginary_Time,stat,Beta,
            #                           numpy.array([ (n+0.5)*Beta/Nmax for n in range(Nmax)]))
            
        GFBloc_ImTime_cython.__init__(self,*self._prepare_init(d))

    #-----------------------------------------------------

    def _plot_(self, OptionsDict):
        """ Plot protocol. OptionsDict can contain : 
             * :param RI: 'R', 'I', 'RI' [ default] 
             * :param x_window: (xmin,xmax) or None [default]
             * :param Name: a string [default ='']. If not '', it remplaces the name of the function just for this plot.
        """
        has_complex_value = False
        M = [x for x in self.mesh]
        return self._plot_base( OptionsDict,  r'$\tau$', lambda name : r'%s$(\tau)$'%name, has_complex_value , M)
 
#-----------------------------------------------------
#  Register the class for HDF_Archive
#-----------------------------------------------------

from pytriqs.Base.Archive.HDF_Archive_Schemes import register_class
register_class (GFBloc_ImTime)


