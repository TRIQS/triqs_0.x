
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

from pytriqs_GF import GF_Statistic, GF_Type, MeshGf
from gf_base import GfBase
import numpy

#-----------------------------------------------------
#  Code Injection
#-----------------------------------------------------

from pytriqs.base.utility.injector import make_injector        # inject new code in the SAME class
from pytriqs_GF import GfLegendre     # the wrapped C++ class.

class __inject (make_injector(GfLegendre), GfBase, GfLegendre):
    """
    A matrix-valued block Green's function described using Legendre coefficients.
    """
    hdf5_scheme_doc = {'data' : "The array of data"}

    def __init__(self, **d):
        """
     The constructor has two variants : you can either provide the mesh
     yourself (see below), or give the parameters to build it.
     All parameters must be given with keyword arguments.

     GfLegendre (indices, beta, statistic, n_legendre_coeffs, data, tail, name, note)

           * ``indices``: a list of indices names of the block
           * ``beta``: the inverse Temperature 
           * ``statistic``: GF_Statistic.Fermion [default] or GF_Statistic.Boson
           * ``n_legendre_coeffs``:  the number of Legendre coefficients to be used
           * ``data``:  a numpy array of dimensions (len(indices),len(indices),n_legendre_coeffs) representing the values of the coefficients.
           * ``tail``:  the tail 
           * ``name``:  a name for the Green's function
           * ``note``:  any string you like...

     If you already have the mesh, you can use a simpler version:

     GfLegendre (indices, mesh, data, tail, name,note)
        
           * ``indices``:  a list of indices names of the block
           * ``mesh``:  a MeshGf object, such that mesh.TypeGF == GF_Type.Imaginary_Legendre
           * ``data``:  a numpy array of dimensions (len(indices),len(indices),n_legendre_coeffs) representing the value of the coefficients.
           * ``tail``:  the tail 
           * ``name``:  a name for the Green's function
           * ``note``:  any string you like...

.. warning::
    The Green function take a **view** of the array data, and a **reference** to the tail.
    """
        if 'mesh' not in d:
            if 'beta' not in d: raise ValueError, "beta not provided"
            beta = float(d['beta'])
            Nmax = d['n_legendre_coeffs'] if 'n_legendre_coeffs' in d else 30
            stat = d['statistic'] if 'statistic' in d else GF_Statistic.Fermion
            sh = 1 if stat == GF_Statistic.Fermion else 0
            d['mesh'] = MeshGf( GF_Type.Imaginary_Legendre, stat, beta,
                                numpy.array(range(Nmax)) )
            for a in ['beta', 'statistic', 'n_legendre_coeffs']:
                if a in d : del d[a]
        else:
            assert d['mesh'].TypeGF==GF_Type.Imaginary_Legendre, "You provided a wrong type of mesh !!"

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

