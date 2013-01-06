
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

__all__ = ['GfReFreq']
from pytriqs_GF import GF_Statistic,GF_Type,TailGf,MeshGf
from gf_base import GfBase
import numpy
from math import pi

#-----------------------------------------------------
#  Code Injection
#-----------------------------------------------------

from pytriqs.base.utility.injector import make_injector        # inject new code in the SAME class
from pytriqs_GF import GfReFreq     # the wrapped C++ class.

class __inject (make_injector(GfReFreq), GfBase, GfReFreq):
    """ 
    A matrix-valued block Green's function in real frequencies.
    """
    def __init__(self,**d):
        """
     The constructor has two variants : you can either provide the mesh in
     real frequencies yourself, or give the parameters to build it.
     All parameters must be given with keyword arguments.

     GfReFreq(indices, beta, statistic, mesh_array, data, tail, name, note)

           * ``indices``:  a list of indices names of the block
           * ``beta``:  Inverse Temperature 
           * ``statistic``:  GF_Statistic.Fermion [default] or GF_Statistic.Boson
           * ``mesh_array``: Grid of frequencies, as a numpy array.
           * ``data``:   A numpy array of dimensions (len(indices),len(indices),len(mesh_array)) representing the value of the Green function on the mesh. 
           * ``tail``:  the tail 
           * ``name``:  a name of the Green's function
           * ``note``:  any string you like...

     If you already have the mesh, you can use a simpler version :

     GfReFreq(indices, mesh, data, tail, name, note)
        
           * ``indices``:  a list of indices names of the block
           * ``mesh``:  a MeshGf object, such that mesh.TypeGF== GF_Type.Real_Frequency 
           * ``data``:   A numpy array of dimensions (len(indices),len(indices),len(mesh_array)) representing the value of the Green function on the mesh.            * ``tail``:  the tail 
           * ``name``:  a name of the Green's function
           * ``note``:  any string you like...

.. warning::
    The Green function take a **view** of the array data, and a **reference** to the tail.

        """
        # construct the mesh if needed
        if 'mesh' not in d : 
            if 'beta' not in d : raise ValueError, "beta not provided"
            if 'mesh_array' not in d : raise ValueError, "mesh_array not provided"
            beta = float(d['beta'])
            stat = d['statistic'] if 'statistic' in d else GF_Statistic.Fermion
            sh = 1 if stat== GF_Statistic.Fermion else 0
            d['mesh'] = MeshGf( GF_Type.Real_Frequency,stat,beta,d['mesh_array'])
            for a in [ 'beta', 'statistic', 'mesh_array'] : 
                if a in d : del d[a]
        else : 
            assert d['mesh'].TypeGF== GF_Type.Real_Frequency, "You provided a wrong type of mesh !!"

        self._init_base__(d)
        self._init_before_injection__(*self._param_for_cons)
        del self._param_for_cons
       
    #-----------------------------------------------------

    def _plot_(self, OptionsDict):
        """ Plot protocol. OptionsDict can contain : 
             * :param RIS: 'R', 'I', 'S', 'RI' [ default] 
             * :param x_window: (xmin,xmax) or None [default]
             * :param name: a string [default ='']. If not '', it remplaces the name of the function just for this plot.
        """
        M = [x for x in self.mesh]
        return self._plot_base( OptionsDict,  r'$\omega$', lambda name : r'%s$(\omega)$'%name, True, M)
 
    #-----------------------------------------------------

    def InverseFourier(self, time_min =None) : 
        """
           Returns a GfReTime containing the Inverse Fourier transform of self
           time_min is the minimal time. By default the time window is centered around 0
        """
        import gf_retime
        (a,b),N = list(self.mesh)[0:2], len(self.mesh)
        om0 = b-a
        if time_min !=None :
            time_max = time_min + 2* pi/om0
        else :
            time_max = pi/om0
            time_min = -time_max
        gt = gf_retime.GfReTime( indices = self.indices, beta = self.beta,
                                           statistic = self.statistic,
                                           time_min = time_min, time_max = time_max, n_time_slices = N )
        gt.setFromInverseFourierOf(self)
        return gt


#-----------------------------------------------------
#  Register the class for HDFArchive
#-----------------------------------------------------

from pytriqs.base.archive.hdf_archive_schemes import register_class
register_class (GfReFreq)
 
 
