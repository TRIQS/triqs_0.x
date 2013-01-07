
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

__all__ = ['GfImTime']
from pytriqs_GF import GF_Statistic,GF_Type,TailGf,MeshGf
from gf_base import GfBase
import numpy
from math import pi

#-----------------------------------------------------
#  Code Injection
#-----------------------------------------------------

from pytriqs.base.utility.injector import make_injector        # inject new code in the SAME class
from pytriqs_GF import GfImTime     # the wrapped C++ class.

class __inject (make_injector(GfImTime), GfBase, GfImTime):
    """ 
    A matrix-valued block Green's function in Matsubara time.
    """
    def __init__(self, **d):
        """
  
     The constructor have two variants : you can either provide the mesh in
     Matsubara frequencies yourself, or give the parameters to build it.
     All parameters must be given with keyword arguments.

     GfImTime(indices, beta, statistic, n_time_points,  data, tail, name, note)
           * ``indices``:  a list of indices names of the block
           * ``beta``:  Inverse Temperature 
           * ``statistic``:  GF_Statistic.Fermion [default] or GF_Statistic.Boson
           * ``n_time_points``  : Number of time slices in the mesh
           * ``data``:   A numpy array of dimensions (len(indices),len(indices),n_time_points) representing the value of the Green function on the mesh. 
           * ``tail``:  the tail 
           * ``name``:  a name of the Green's function
           * ``note``:  any string you like...

     If you already have the mesh, you can use a simpler version :

     GfImTime (indices, mesh, data, tail, name,note)
        
           * ``indices``:  a list of indices names of the block
           * ``mesh``:  a MeshGf object, such that Mesh.TypeGF== GF_Type.Imaginary_Time 
           * ``data``:   A numpy array of dimensions (len(indices),len(indices),n_time_points) representing the value of the Green function on the mesh. 
           * ``tail``:  the tail 
           * ``name``:  a name of the Green's function
           * ``note``:  any string you like...

 .. warning::
    The Green function take a **view** of the array data, and a **reference** to the tail.

         """
        # construct the mesh if needed
        if 'mesh' not in d : 
            if 'beta' not in d : raise ValueError, "beta not provided"
            beta = float(d['beta'])
            Nmax = d['n_time_points'] if 'n_time_points' in d else 10000
            stat = d['statistic'] if 'statistic' in d else GF_Statistic.Fermion
            sh = 1 if stat== GF_Statistic.Fermion else 0
            d['mesh'] = MeshGf( GF_Type.Imaginary_Time,stat,beta,
                                       numpy.array([ (n+0.5)*beta/Nmax for n in range(Nmax)]))
            for a in [ 'beta', 'statistic', 'n_time_points'] : 
                if a in d : del d[a]
        else : 
            assert d['mesh'].TypeGF==GF_Type.Imaginary_Time, "You provided a wrong type of mesh !!"
            
        self._init_base__(d)
        self._init_before_injection__(*self._param_for_cons)
        del self._param_for_cons
                
    #-----------------------------------------------------

    def _plot_(self, OptionsDict):
        """ Plot protocol. OptionsDict can contain : 
             * :param RI: 'R', 'I', 'RI' [ default] 
             * :param x_window: (xmin,xmax) or None [default]
             * :param name: a string [default ='']. If not '', it remplaces the name of the function just for this plot.
        """
        has_complex_value = False
        M = [x for x in self.mesh]
        return self._plot_base( OptionsDict,  r'$\tau$', lambda name : r'%s$(\tau)$'%name, has_complex_value , M)
 
#-----------------------------------------------------
#  Register the class for HDFArchive
#-----------------------------------------------------

from pytriqs.base.archive.hdf_archive_schemes import register_class
register_class (GfImTime)


