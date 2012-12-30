################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by L. Boehnke, M. Ferrero, O. Parcollet
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

__all__ = ['GFBloc_ImLegendre']
import numpy
from gf import GFBloc_ReTime_cython, MeshReTime 
from GFBloc_general import _GFBloc_general 

class GFBloc_ImLegendre (GFBloc_ImLegendre_cython, _GFBloc_general):
    """
    A matrix-valued block Green's function described using Legendre coefficients.
    """
    hdf5_scheme_doc = {'data' : "The array of data"}

    def __init__(self, **d):
        """
     The constructor has two variants : you can either provide the mesh
     yourself (see below), or give the parameters to build it.
     All parameters must be given with keyword arguments.

     GFBloc_ImLegendre(Indices, Beta, Statistic, NLegendreCoeffs, Data, Tail, Name, Note)

           * ``Indices``: a list of indices names of the block
           * ``Beta``: the inverse Temperature 
           * ``Statistic``: GF_Statistic.Fermion [default] or GF_Statistic.Boson
           * ``NLegendreCoeffs``:  the number of Legendre coefficients to be used
           * ``Data``:  a numpy array of dimensions (len(Indices),len(Indices),NLegendreCoeffs) representing the values of the coefficients.
           * ``Tail``:  the tail 
           * ``Name``:  a name for the Green's function
           * ``Note``:  any string you like...

     If you already have the mesh, you can use a simpler version:

     GFBloc_ImLegendre(Indices, Mesh, Data, Tail, Name,Note)
        
           * ``Indices``:  a list of indices names of the block
           * ``Mesh``:  a MeshGF object, such that Mesh.TypeGF == GF_Type.Imaginary_Legendre
           * ``Data``:  a numpy array of dimensions (len(Indices),len(Indices),NLegendreCoeffs) representing the value of the coefficients.
           * ``Tail``:  the tail 
           * ``Name``:  a name for the Green's function
           * ``Note``:  any string you like...

.. warning::
    The Green function take a **view** of the array Data, and a **reference** to the Tail.
    """
       # construct the mesh if needed
        if 'Mesh' not in d : 
            if 'Beta' not in d : raise ValueError, "Beta not provided"
            Beta = float(d.pop('Beta'))
            Nmax = d.pop('NLegendreCoeffs',30)
            stat = d.pop('Statistic','F') # GF_Statistic.Fermion
            sh = 1 if stat== 'F' else 0 # GF_Statistic.Fermion else 0
            d['Mesh'] = MeshLegendre(Beta,'F',Nmax)
            #d['Mesh'] = MeshGF( GF_Type.Imaginary_Legendre, stat, Beta, numpy.array(range(Nmax)) )

        GFBloc_ImLegendre_cython.__init__(self,*self._prepare_init(d))
                
    #-----------------------------------------------------

    def _plot_(self, OptionsDict):
        """ Plot protocol. OptionsDict can contain : 
             * :param RIS: 'R', 'I', 'S', 'RI' [ default] 
             * :param x_window: (xmin,xmax) or None [default]
             * :param Name: a string [default ='']. If not '', it remplaces the name of the function just for this plot.
        """
        M = [x for x in self.mesh]
        if "RI" not in OptionsDict: OptionsDict["RI"] = "R"
        return self._plot_base( OptionsDict, r'$l$', lambda name : r'%s$_l$'%name, True, M)

    #-----------------------------------------------------

    def copy_and_truncate(self, Nleg):
        """ Copies the Legendre Green's function and truncates it down

        Parameters
        ----------
        Nleg : int
          remaining number of Legendre coefficients after truncation
        """
        new_g = self.__class__(IndicesL = self._IndicesL,
                               IndicesR = self._IndicesR,
                               Beta = self.Beta,
                               Statistic = self.Statistic,
                               NLegendreCoeffs = Nleg,
                               Name = self.Name, Note = self.Note)
        new_g._data.array[:,:,:] = self._data.array[:,:,0:Nleg]
        new_g.determine_tail()
        return new_g


#-----------------------------------------------------
#  Register the class for HDF_Archive
#-----------------------------------------------------

from pytriqs.base.Archive.HDF_Archive_Schemes import register_class
register_class (GFBloc_ImLegendre)


