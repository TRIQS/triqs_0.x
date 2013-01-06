
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
import numpy
from math import pi
from gf import GfReFreq_cython, MeshReFreq  
from GFBloc_general import _GFBloc_general 

class GfReFreq (GfReFreq_cython, _GFBloc_general):
    """ 
    A matrix-valued block Green's function in real frequencies.
    """
    def __init__(self,**d):
        """
     The constructor has two variants : you can either provide the mesh in
     real frequencies yourself, or give the parameters to build it.
     All parameters must be given with keyword arguments.

     GfReFreq(Indices, Beta, Statistic, MeshArray, Data, Tail, Name,Note)

           * ``Indices``:  a list of indices names of the block
           * ``Beta``:  Inverse Temperature 
           * ``Statistic``:  GF_Statistic.Fermion [default] or GF_Statistic.Boson
           * ``MeshArray``: Grid of frequencies, as a numpy array.
           * ``Data``:   A numpy array of dimensions (len(Indices),len(Indices),len(MeshArray)) representing the value of the Green function on the mesh. 
           * ``Tail``:  the tail 
           * ``Name``:  a name of the GF
           * ``Note``:  any string you like...

     If you already have the mesh, you can use a simpler version :

     GfReFreq(Indices, Mesh, Data, Tail, Name,Note)
        
           * ``Indices``:  a list of indices names of the block
           * ``Mesh``:  a MeshGf object, such that Mesh.TypeGF== GF_Type.Real_Frequency 
           * ``Data``:   A numpy array of dimensions (len(Indices),len(Indices),len(MeshArray)) representing the value of the Green function on the mesh.            * ``Tail``:  the tail 
           * ``Name``:  a name of the GF
           * ``Note``:  any string you like...

.. warning::
    The Green function take a **view** of the array Data, and a **reference** to the Tail.

        """
        # construct the mesh if needed
        if 'Mesh' not in d : 
            if 'Beta' not in d : raise ValueError, "Beta not provided"
            Beta = float(d.pop('Beta'))
            Nmax = d.pop('Nmax',1025)
            stat = d.pop('Statistic','F') # GF_Statistic.Fermion
            sh = 1 if stat== 'F' else 0 # GF_Statistic.Fermion else 0
            d['Mesh'] = MeshRealFrequency(Beta,Nmax)
            #if 'MeshArray' not in d : raise ValueError, "MeshArray not provided"
            #d['Mesh'] = MeshGf( GF_Type.Real_Frequency,stat,Beta,d['MeshArray'])

        GfReFreq_cython.__init__(self,*self._prepare_init(d))
  
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
           Returns a GfReTime containing the Inverse Fourier transform of self
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
        gt = gf_retime.GfReTime( indices = self.indices,beta = self.beta,
                                          Statistic = self.Statistic,
                                          TimeMin = TimeMin, TimeMax = TimeMax,NTimeSlices = N )
        gt.setFromInverseFourierOf(self)
        return gt


#-----------------------------------------------------
#  Register the class for HDFArchive
#-----------------------------------------------------

from pytriqs.base.archive.hdf_archive_schemes import register_class
register_class (GfReFreq)
 
 
