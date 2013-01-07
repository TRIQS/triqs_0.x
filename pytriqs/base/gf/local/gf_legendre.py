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

__all__ = ['GfLegendre']
import numpy
from gf import GfReTime_cython, MeshReTime 
from GFBloc_general import _GFBloc_general 

class GfLegendre (GfLegendre_cython, _GFBloc_general):
    """
    A matrix-valued block Green's function described using Legendre coefficients.
    """
    hdf5_scheme_doc = {'data' : "The array of data"}

    def __init__(self, **d):
        """
     The constructor has two variants : you can either provide the mesh
     yourself (see below), or give the parameters to build it.
     All parameters must be given with keyword arguments.

     GfLegendre(indices, beta, statistic, n_legendre_coeffs, data, tail, name, note)

           * ``indices``: a list of indices names of the block
           * ``beta``: the inverse Temperature 
           * ``statistic``: 'F' or 'B'
           * ``n_legendre_coeffs``:  the number of Legendre coefficients to be used
           * ``data``:  a numpy array of dimensions (len(indices),len(indices),n_legendre_coeffs) representing the values of the coefficients.
           * ``tail``:  the tail 
           * ``name``:  a name for the Green's function
           * ``note``:  any string you like...

     If you already have the mesh, you can use a simpler version:

     GfLegendre(indices, mesh, data, tail, name, note)
        
           * ``indices``:  a list of indices names of the block
           * ``mesh``:  a MeshGf object, such that mesh.TypeGF == GF_Type.Imaginary_Legendre
           * ``data``:  a numpy array of dimensions (len(indices),len(indices),n_legendre_coeffs) representing the value of the coefficients.
           * ``tail``:  the tail 
           * ``name``:  a name for the Green's function
           * ``note``:  any string you like...

.. warning::
    The Green function take a **view** of the array data, and a **reference** to the tail.
    """
       # construct the mesh if needed
        if 'mesh' not in d : 
            if 'beta' not in d : raise ValueError, "beta not provided"
            beta = float(d.pop('beta'))
            n_max = d.pop('n_legendre_coeffs',30)
            stat = d.pop('statistic','F') # GF_statistic.Fermion
            sh = 1 if stat== 'F' else 0 # GF_statistic.Fermion else 0
            d['mesh'] = MeshLegendre(beta,'F',n_max)
            #d['mesh'] = MeshGf( GF_Type.Imaginary_Legendre, stat, beta, numpy.array(range(n_max)) )

        GfLegendre_cython.__init__(self,*self._prepare_init(d))
                
    #-----------------------------------------------------

    def _plot_(self, OptionsDict):
        """ Plot protocol. OptionsDict can contain : 
             * :param RIS: 'R', 'I', 'S', 'RI' [ default] 
             * :param x_window: (xmin,xmax) or None [default]
             * :param name: a string [default ='']. If not '', it remplaces the name of the function just for this plot.
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
        new_g = self.__class__(indicesL = self._indicesL,
                               indicesR = self._indicesR,
                               beta = self.beta,
                               statistic = self.statistic,
                               n_legendre_coeffs = Nleg,
                               name = self.name, note = self.note)
        new_g._data.array[:,:,:] = self._data.array[:,:,0:Nleg]
        new_g.determine_tail()
        return new_g


#-----------------------------------------------------
#  Register the class for HDFArchive
#-----------------------------------------------------

from pytriqs.base.archive.hdf_archive_schemes import register_class
register_class (GfLegendre)


