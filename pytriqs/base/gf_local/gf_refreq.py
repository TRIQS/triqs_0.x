
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

__all__ = ['GFBloc_ReFreq']
from pytriqs_GF import GF_Statistic,GF_Type,TailGF,MeshGF
from gf_base import gf_base
from gf_concept import gf_concept
import numpy
from math import pi

#-----------------------------------------------------
#  Code Injection
#-----------------------------------------------------

from pytriqs.base.utility.injector import make_injector        # inject new code in the SAME class
from pytriqs_GF import GFBloc_ReFreq     # the wrapped C++ class.

class __inject (make_injector(GFBloc_ReFreq) ,gf_concept,gf_base, GFBloc_ReFreq):
    """ 
    A matrix-valued block Green's function in real frequencies.
    """
    def __init__(self,**d):
        """
     The constructor has two variants : you can either provide the mesh in
     real frequencies yourself, or give the parameters to build it.
     All parameters must be given with keyword arguments.

     GFBloc_ReFreq(Indices, Beta, Statistic, MeshArray, Data, Tail, Name,Note)

           * ``Indices``:  a list of indices names of the block
           * ``Beta``:  Inverse Temperature 
           * ``Statistic``:  GF_Statistic.Fermion [default] or GF_Statistic.Boson
           * ``MeshArray``: Grid of frequencies, as a numpy array.
           * ``Data``:   A numpy array of dimensions (len(Indices),len(Indices),len(MeshArray)) representing the value of the Green function on the mesh. 
           * ``Tail``:  the tail 
           * ``Name``:  a name of the GF
           * ``Note``:  any string you like...

     If you already have the mesh, you can use a simpler version :

     GFBloc_ReFreq(Indices, Mesh, Data, Tail, Name,Note)
        
           * ``Indices``:  a list of indices names of the block
           * ``Mesh``:  a MeshGF object, such that Mesh.TypeGF== GF_Type.Real_Frequency 
           * ``Data``:   A numpy array of dimensions (len(Indices),len(Indices),len(MeshArray)) representing the value of the Green function on the mesh.            * ``Tail``:  the tail 
           * ``Name``:  a name of the GF
           * ``Note``:  any string you like...

.. warning::
    The Green function take a **view** of the array Data, and a **reference** to the Tail.

        """
        # construct the mesh if needed
        if 'Mesh' not in d : 
            if 'Beta' not in d : raise ValueError, "Beta not provided"
            if 'MeshArray' not in d : raise ValueError, "MeshArray not provided"
            Beta = float(d['Beta'])
            stat = d['Statistic'] if 'Statistic' in d else GF_Statistic.Fermion
            sh = 1 if stat== GF_Statistic.Fermion else 0
            d['Mesh'] = MeshGF( GF_Type.Real_Frequency,stat,Beta,d['MeshArray'])
            for a in [ 'Beta', 'Statistic', 'MeshArray'] : 
                if a in d : del d[a]
        else : 
            assert d['Mesh'].TypeGF== GF_Type.Real_Frequency, "You provided a wrong type of mesh !!"

        self._init_base__(d)
        self._init_before_injection__(*self._param_for_cons)
        del self._param_for_cons
       
    #-----------------------------------------------------

    def _plot_(self, OptionsDict):
        """ Plot protocol. OptionsDict can contain : 
             * :param RIS: 'R', 'I', 'S', 'RI' [ default] 
             * :param x_window: (xmin,xmax) or None [default]
             * :param Name: a string [default ='']. If not '', it remplaces the name of the function just for this plot.
        """
        M = [x for x in self.mesh]
        return self._plot_base( OptionsDict,  r'$\omega$', lambda name : r'%s$(\omega)$'%name, True, M)
 
    #-----------------------------------------------------

    def InverseFourier(self, TimeMin =None) : 
        """
           Returns a GFBloc_ReTime containing the Inverse Fourier transform of self
           TimeMin is the minimal time. By default the time window is centered around 0
        """
        import gf_retime
        (a,b),N = list(self.mesh)[0:2], len(self.mesh)
        om0 = b-a
        if TimeMin !=None :
            TimeMax = TimeMin + 2* pi/om0
        else :
            TimeMax = pi/om0
            TimeMin = -TimeMax
        gt = gf_retime.GFBloc_ReTime( Indices = self.Indices,Beta = self.Beta,
                                          Statistic = self.Statistic,
                                          TimeMin = TimeMin, TimeMax = TimeMax,NTimeSlices = N )
        gt.setFromInverseFourierOf(self)
        return gt


#-----------------------------------------------------
#  Register the class for HDFArchive
#-----------------------------------------------------

from pytriqs.base.archive.hdf_archive_schemes import register_class
register_class (GFBloc_ReFreq)
 
 
