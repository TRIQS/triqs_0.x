
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

from scipy.optimize import leastsq
import numpy as np

class general_fit : 

    def __init__ (self, X,Y, f, f_str = None ) : 
        self.f, self.f_str = f,f_str
        errfunc = lambda x :  np.abs ( f(X,*x)  - Y[:])
        self.param, success = leastsq(errfunc, (1,1) )

    def __str__ (self) : return self.f_str(*self.param) if self.f_str else 'Fit'

    def __call__ (self,omega) : return self.f(omega, *self.param).imag



linear_fit = lambda X,Y :general_fit (X,Y, 
                                      f = lambda Omega, a, b : a * Omega *1j  + b,
                                      f_str = lambda a,b : r"$%s *\omega + %s$"%(a,b) )

onelevel_fit = lambda X,Y :general_fit (X,Y,
                                        f= lambda Omega, a, b : 1/(a * Omega * 1j + b), 
                                        f_str = lambda a,b : r"${1}/(%s *\omega + %s)$"%(a,b)   )

if __name__ == "__main__" : 

    from pytriqs.Base.GF_Local import *
    from pytriqs.Base.GF_Local.Descriptors import Omega
    g = GFBloc_ImFreq(Indices = [1], Beta = 50, NFreqMatsubara = 1000, Name = "g")
    g <<= inverse( Omega + 0.5 )

    X,Y = g.x_data_view (x_window = (0,0.2), flatten_y = True ) 
    fit = linear_fit    ( X,Y )
    fit2 = onelevel_fit ( X,Y )
    
    from pytriqs.Base.Plot.MatplotlibInterface import plot
    
    plot (g, '-o', x_window = (0,10) ) 
    plot (fit , '-x', x_window = (0,0.5) )
    plot (fit2 , '-x', x_window = (0,5) )
    
