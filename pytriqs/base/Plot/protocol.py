
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

import numpy

def clip_array (X, xmin, xmax) : 
    """
     Given : 
      - X a 1d numpy array of shape (L), of ordered values
        or in fact any generator of ordered values.
      - xmin, xmax
      it returns the slice sl such that 
        * X[sl] is in [xmin, xmax]
    """
    try : 
        low  = (i for i,x in enumerate(X) if not( x < xmin)  ).next() 
    except StopIteration :
        low = 0

    try :
        high = (i for i,x in enumerate(X) if x > xmax ).next()
        r = slice(low,high)
    except StopIteration :
        r = slice(low,)

    return r

def plot_protocol_apply(ob, OptionsDict, xlims) : 
    """
    Given an object ob that supports the plot protocol, it applies the protocol
    or emulate it
    """
 
    if hasattr(ob, '_plot_') : return ob._plot_(OptionsDict)
    elif callable(ob) :
        NPoints = OptionsDict.pop('NPoints',100)        
        rx = OptionsDict.pop('x_window',None ) # consume it
        xmin,xmax = rx if rx else xlims()
        X = numpy.arange(xmin,xmax, (xmax- xmin)/float(NPoints))
        Y = numpy.array( [ob(x) for x in X] )
    else :
        try : # generator x,y
            X, Y = izip (*ob)
        except : 
            raise RuntimeError, "Object can not be plotted"

    Name = OptionsDict.pop('Name','') 
    if not Name : Name = str(ob)
    if numpy.iscomplexobj(Y) : 
        return( [  {'Type' : "XY", 'xdata':X, 'ydata':Y.real, 'label': "Re " + Name} , 
                   {'Type' : "XY", 'xdata':X, 'ydata':Y.imag, 'label': "Im " + Name} ] )
    else:
        return( [ {'Type' : "XY", 'xdata':X, 'ydata':Y, 'label': Name} ] )
 
