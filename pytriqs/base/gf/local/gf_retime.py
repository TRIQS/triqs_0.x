
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

__all__ = ['GfReTime']
import numpy
from math import pi
from gf import GfReTime_cython, MeshReTime 
from GFBloc_general import _GFBloc_general 

class GfReTime (GfReTime_cython, _GFBloc_general):
    """ 
    A matrix-valued block Green's function in real time.
    """
    def __init__(self, **d):
        """
    The constructor has two variants : you can either provide the mesh in
    real time yourself, or give the parameters to build it.
    All parameters must be given with keyword arguments.

    GfReTime (Indices, Beta, Statistic, NTimeSlices, TimeMin, TimeMax,  Data, Tail, Name,Note)

           * ``Indices``:  a list of indices names of the block
           * ``Beta``:  Inverse Temperature 
           * ``Statistic``:  GF_Statistic.Fermion [default] or GF_Statistic.Boson
           * ``NTimeSlices``  : Number of time slices
           * ``TimeMin,TimeMax``  : The time window
           * ``Data``:   A numpy array of dimensions (len(Indices),len(Indices),NTimeSlices) representing the value of the Green function on the mesh. 
           * ``Tail``:  the tail 
           * ``Name``:  a name of the GF
           * ``Note``:  any string you like...

    If you already have the mesh, you can use a simpler version :

    GfReTime (Indices, Mesh, Data, Tail, Name,Note)
        
           * ``Indices``:  a list of indices names of the block
           * ``Mesh``:  a MeshGf object, such that Mesh.TypeGF== GF_Type.Real_Time 
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
            Nmax = d.pop('NTimeSlices',1024)
            assert Nmax%2 ==0, "Use an even number of slices"
            timeMin = d.pop('TimeMin',-10)
            timeMax = d.pop('TimeMax', 10)
            dt = float (timeMax - timeMin)/ Nmax
            d['Mesh'] = MeshRealTime( GF_Type.Real_Time,stat,Beta,
                                numpy.array([ timeMin + (n+0.5)*dt for n in range(Nmax)]))

        self.TimeMin, self.TimeMax, self.Npts = min(d['Mesh']), max(d['Mesh']), len(d['Mesh'])
        dt = (self.TimeMax - self.TimeMin)/(self.Npts -1)
        self.TimeMin -= dt/2 ;self.TimeMax += dt/2 
        # ????? duplicated in C++... not very clean, but ok.

        GfReTime_cython.__init__(self,*self._prepare_init(d))
        
    #-----------------------------------------------------
    
    def Fourier(self):
        """Returns a GfReFreq containing the Fourier transform of self"""
        import gf_refreq
        om0 = 2*pi/(self.TimeMax - self.TimeMin)
        N = self.Npts
        gw = gf_refreq.GfReFreq(indices = self.indices,beta = self.beta,
                                         Statistic = self.Statistic,
                                         MeshArray = numpy.array([ om0*i for i in range (- (N/2),N/2)]))
        gw.setFromFourierOf(self)
        return gw
                
    #-----------------------------------------------------

    def _plot_(self, OptionsDict):
        """ Plot protocol. OptionsDict can contain : 
             * :param RI: 'R', 'I', 'RI' [ default] 
             * :param x_window: (xmin,xmax) or None [default]
             * :param Name: a string [default ='']. If not '', it remplaces the name of the function just for this plot.
        """
        M = [x for x in self.mesh]
        return self._plot_base( OptionsDict,  r'$t$', lambda name : r'%s$(t)$'%name, True, M)
 
#-----------------------------------------------------
#  Register the class for HDFArchive
#-----------------------------------------------------

from pytriqs.base.archive.hdf_archive_schemes import register_class
register_class (GfReTime)



 
