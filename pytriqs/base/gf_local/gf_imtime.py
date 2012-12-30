
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
from pytriqs_GF import GF_Statistic,GF_Type,TailGF,MeshGF
from gf_base import gf_base
from gf_concept import gf_concept
import numpy
from math import pi

#-----------------------------------------------------
#  Code Injection
#-----------------------------------------------------

from pytriqs.base.utility.Injector import make_injector        # inject new code in the SAME class
from pytriqs_GF import GFBloc_ImTime     # the wrapped C++ class.

class __inject (make_injector(GFBloc_ImTime) ,gf_concept, gf_base, GFBloc_ImTime):
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
            Beta = float(d['Beta'])
            Nmax = d['NTimeSlices'] if 'NTimeSlices' in d else 10000
            stat = d['Statistic'] if 'Statistic' in d else GF_Statistic.Fermion
            sh = 1 if stat== GF_Statistic.Fermion else 0
            d['Mesh'] = MeshGF( GF_Type.Imaginary_Time,stat,Beta,
                                       numpy.array([ (n+0.5)*Beta/Nmax for n in range(Nmax)]))
            for a in [ 'Beta', 'Statistic', 'NTimeSlices'] : 
                if a in d : del d[a]
        else : 
            assert d['Mesh'].TypeGF==GF_Type.Imaginary_Time, "You provided a wrong type of mesh !!"
            
        self._init_base__(d)
        self._init_before_injection__(*self._param_for_cons)
        del self._param_for_cons
                
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

from pytriqs.base.archive.hdf_archive_schemes import register_class
register_class (GFBloc_ImTime)


