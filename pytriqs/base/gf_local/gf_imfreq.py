
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

__all__ = ['GfImFreq']
from pytriqs_GF import GF_Statistic,GF_Type,TailGf,MeshGf
from gf_base import GfBase
import numpy
from math import pi

#-----------------------------------------------------
#  Code Injection
#-----------------------------------------------------

from pytriqs.base.utility.injector import make_injector        # inject new code in the SAME class
from pytriqs_GF import GfImFreq     # the wrapped C++ class.

class __inject (make_injector(GfImFreq), GfBase, GfImFreq):
    """ 
    A matrix-valued block Green's function in Matsubara frequencies.
    """
    hdf5_scheme_doc = {'data' : "The array of data"}
    
    def __init__(self, **d):
        """
     The constructor have two variants : you can either provide the mesh in
     Matsubara frequencies yourself, or give the parameters to build it.
     All parameters must be given with keyword arguments.

     GfImFreq(indices, beta, statistic, n_matsubara, data, tail, name, note)

           * ``indices``:  a list of indices names of the block
           * ``beta``:  Inverse Temperature 
           * ``statistic``:  GF_Statistic.Fermion [default] or GF_Statistic.Boson
           * ``n_matsubara``:  Number of Matsubara frequencies
           * ``data``:   A numpy array of dimensions (len(indices),len(indices),n_matsubara) representing the value of the Green function on the mesh. 
           * ``tail``:  the tail 
           * ``name``:  a name of the Green's function
           * ``note``:  any string you like...

     If you already have the mesh, you can use a simpler version :

     GfImFreq(indices, mesh, data, tail, name, note)
        
           * ``indices``:  a list of indices names of the block
           * ``mesh``:  a MeshGf object, such that mesh.TypeGF== GF_Type.Imaginary_Frequency 
           * ``data``:   A numpy array of dimensions (len(indices),len(indices),n_matsubara) representing the value of the Green function on the mesh. 
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
            Nmax = d['n_matsubara'] if 'n_matsubara' in d else 1025
            stat = d['statistic'] if 'statistic' in d else GF_Statistic.Fermion
            sh = 1 if stat== GF_Statistic.Fermion else 0
            d['mesh'] = MeshGf( GF_Type.Imaginary_Frequency,stat,beta,
                           numpy.array([ (2*n+sh)*pi/beta for n in range(Nmax)]))
            for a in [ 'beta', 'statistic', 'n_matsubara'] : 
                if a in d : del d[a]
        else : 
            assert d['mesh'].TypeGF== GF_Type.Imaginary_Frequency, "You provided a wrong type of mesh !!"

        self._init_base__(d)
        self._init_before_injection__(*self._param_for_cons)
        del self._param_for_cons

    #-----------------------------------------------------

    def _plot_(self, opt_dict):
        """ Plot protocol. opt_dict can contain : 
             * :param RIS: 'R', 'I', 'S', 'RI' [ default] 
             * :param x_window: (xmin,xmax) or None [default]
             * :param name: a string [default ='']. If not '', it remplaces the name of the function just for this plot.
        """
        M = [x.imag for x in self.mesh]
        return self._plot_base(opt_dict,  r'$\omega_n$', lambda name : r'%s$(i\omega_n)$'%name, True, M)
    
    #-----------------------------------------------------

    def fit_tail(self, fixed_coef, order_max, fit_start, fit_stop, replace_tail = True):
       """
       Fit the tail of the Green's function
       Input:
         - fixed_coef: a 3-dim array of known coefficients for the fit starting from the order -1
         - order_max: highest order in the fit
         - fit_start, fit_stop: fit the data between fit_start and fit_stop
       Output:
         On output all the data above fit_start is replaced by the fitted tail
         and the new moments are included in the Green's function
       """
       from scipy.optimize import leastsq

       # Turn known_coefs into a numpy array if ever it is not already the case
       #known_coef = numpy.array(fixed_coef)
       known_coef = fixed_coef

       # Check the shape of the known_coef
       #assert(known_coef.shape[0:2] == (self.N1,self.N2))

       # Change the OrderMax
       # It is assumed that any known_coef will start at order -1
       self._tail.changeOrderMax(order_max)

       # Fill up two arrays with the frequencies and values over the range of interest
       ninit, nstop = 0, -1
       x = []
       for om in self.mesh:
         if (om.imag < fit_start): ninit = ninit+1
         if (om.imag <= fit_stop): nstop = nstop+1
         if (om.imag <= fit_stop and om.imag >= fit_start): x += [om]
       omegas = numpy.array(x)
       #values = self._data[:,:,ninit:nstop+1].array
       values = self._data.array[:,:,ninit:nstop+1]

       # Loop over the indices of the Green's function
       for n1,indR in enumerate(self._IndicesR):
         for n2,indL in enumerate(self._IndicesL):
           
           # Construct the part of the fitting functions which is known
           f_known = numpy.zeros((len(omegas)),numpy.complex)
           for order in range(len(known_coef[n1][n2])):
             f_known += known_coef[n1][n2][order]*omegas**(1-order)
	   
           # Compute how many free parameters we have and give an initial guess
           len_param = order_max-len(known_coef[n1][n2])+2
           p0 = len_param*[1.0]
	   
           # This is the function to be minimized, the diff between the original
           # data in values and the fitting function
           def fct(p):
             y_fct = 1.0*f_known
             for order in range(len_param):
               y_fct += p[order]*omegas**(1-len(known_coef[n1][n2])-order)
             y_fct -= values[n1,n2,:]
             return abs(y_fct)

           # Now call the minimizing function
           sol = leastsq(fct, p0, maxfev=1000*len_param)

           # Put the known and the new found moments in the tail
           for order in range(len(known_coef[n1][n2])):
             self._tail[order-1][indR,indL] = [[ known_coef[n1][n2][order] ]]
           for order, moment in enumerate(sol[0]):
             self._tail[len(known_coef[n1][n2])+order-1][indR,indL] = [[ moment ]]

       # Replace then end of the Green's function by the tail
       if replace_tail: self.replace_by_tail(ninit);


#-----------------------------------------------------
#  Register the class for HDFArchive
#-----------------------------------------------------

from pytriqs.base.archive.hdf_archive_schemes import register_class
register_class (GfImFreq)


