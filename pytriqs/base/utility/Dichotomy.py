
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

import pytriqs.base.utility.MPI as MPI

def Dichotomy(Function, xinit,yvalue,Precision_on_y,Delta_x, MaxNbreLoop=1000, xname="", yname="",verbosity=1):
    """
    Solver Function(x) = yvalue.
    
    Arguments :
      - Function : function (real valued) to be solved by dichotomy
      - xinit : Init value for x. On success, returns the new value of x
      - yvalue : 
      - precision : calculation stops for abs(f(x) - yvalue)<precision
      - MaxNbreLoop : maximum number of loops before failure. Default is 1000
      - xname, yname : name of the variable x, y for the report
      - verbosity : verbosity level.

    Returns :
       - A tuple (x,y). x is the value found, y is f(x).
       - (None,None) if the calculation failed.

    """
    def sign(x):
        if x>0.0 : return 1
        if x<0.0 : return -1
        return 0
    
    MPI.report("Dichotomy adjustment of %(xname)s to obtain %(yname)s = %(yvalue)f +/- %(Precision_on_y)f"%locals() )
    PR = "    "
    if xname=="" or yname==""  : verbosity = max(verbosity,1)
    x=xinit;Delta_x= abs(Delta_x)

    # First find the bounds
    y1 = Function(x)
    eps = sign(y1-yvalue)
    x1=x;y2=y1;x2=x1
    nbre_loop=0
    while (nbre_loop<= MaxNbreLoop) and (y2-yvalue)*eps>0 and abs(y2-yvalue)>Precision_on_y :
        nbre_loop +=1
        x2 -=  eps*Delta_x
        y2 = Function(x2)
        if xname!="" and verbosity>2:
            MPI.report("%(PR)s%(xname)s = %(x2)f  \n%(PR)s%(yname)s = %(y2)f"%locals())

    MPI.report("%(PR)s%(x1)f < %(xname)s < %(x2)f"%locals())
    MPI.report("%(PR)s%(y1)f < %(yname)s < %(y2)f"%locals())

    # Now mu is between mu1 and mu2
    yfound = y2
    # We found bounds. What if the next loop is never run ?
    # i.e. x1 or x2 are close to the solution
    # we have to know which one is the best .... 
    if abs(y1-yvalue)< abs(y2-yvalue) :
        x=x1
    else:
        x=x2
        
    #Now let's refine our mu....
    while (nbre_loop<= MaxNbreLoop) and (abs(yfound-yvalue)>Precision_on_y) :
        nbre_loop +=1
        x = x1  + (x2 - x1) * (yvalue - y1)/(y2-y1)
        yfound = Function(x)
        if (y1-yvalue)*(yfound - yvalue)>0 : 
            x1 = x; y1=yfound
        else :
            x2= x;y2=yfound;
        if verbosity>2 :
            MPI.report("%(PR)s%(x1)f < %(xname)s < %(x2)f"%locals())
            MPI.report("%(PR)s%(y1)f < %(yname)s < %(y2)f"%locals())
    if abs(yfound - yvalue) < Precision_on_y :
        if verbosity>0:
            MPI.report("%(PR)s%(xname)s found in %(nbre_loop)d iterations : "%locals())
            MPI.report("%(PR)s%(yname)s = %(yfound)f;%(xname)s = %(x)f"%locals())
        return (x,yfound)
    else : 
        if verbosity>0:
            MPI.report("%(PR)sFAILURE to adjust %(xname)s  to the value %(yvalue)f after %(nbre_loop)d iterations."%locals())
        return (None,None)
    
