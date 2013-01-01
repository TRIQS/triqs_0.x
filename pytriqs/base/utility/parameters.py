
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

from types import *
import pytriqs.base.utility.mpi as mpi

class Parameters(dict) : 
    """
    Storage of parameters for a script, a particular class, etc...
    Parameters is a dict with 2 modifications :    
      - It is constructed from various sources : a dict, a file of various formats.
      - It has a method "check" to check and complete the parameters with 
        required, optional parameters, check contraints.

    mpi : By default, the contruction is made on the master, and then bcasted.
          See constructor. 
    """
    def __init__(self, S, BuildMasterOnly =True) :
        """
        Construction : with named arguments only.
        Possible constructors from source S : 
           - Parameters(d) : d is a dict to be copied (no deepcopy, just updated).
           - Parameters(f) : f is a string containing the path of a file
             The type of the file is determined from the extension: 
                - .py, .txt : the file is executed in the Parameters
                - .xml : to be written
        mpi : if BuildMasterOnly is true, the contruction is only made on the master, 
              the result is then bcasted to all the nodes.
              Otherwise, it is done on all nodes (not recommended to read files).
        """
        if mpi.IS_MASTER_NODE() or not BuildMasterOnly : 
            if type(S) == type(''):
               # detect the type of the file
               try : 
                  extension = S.split('.')[-1].lower()
               except IndexError: 
                  raise ValueError, "I am lost : I can not determine the extension of the file !"
               if extension in ['py','txt'] : execfile(S,{},self)
               else : raise ValueError, "Extension of the file not recognized"
            else : # S is therefore a dict 
                try : 
                  self.update(S)
                except : 
                  print "Error in Parameter constructor. Is the source an iterable ?"
                  raise
        # end of master only
        if BuildMasterOnly : mpi.bcast(self) # bcast it on the nodes

# The next functions are not in methods, since we may want to use them 
# on true dictionnaries...

#----------------------------------------------------------------

def parameters_not_in_union_of_dicts (DIC, *d_list):
    """ Given a list of dictionnaries d_list, it checks whether all
    parameters' names are in the union of the keys of those dicts.
    If not, it return the *set* of the extra parameters (which are not in
    this union)"""

    return set(DIC.keys()) - reduce (lambda x,y :  x | set(y.keys()) , d_list, set([]) )
#----------------------------------------------------------------

def check_no_parameters_not_in_union_of_dicts (DIC, *d_list):
    """ Checks that parameters_not_in_union_of_dicts returns void or report and
    raise"""

    p = parameters_not_in_union_of_dicts(DIC,*d_list)
    if p : 
        print "Error. Some parameters are superfluous"
        for x in p : 
          print " *  %s"%x
        raise ValueError, ""

#----------------------------------------------------------------
 
def check (DIC, Required, Optional, Deprecated = (), Constraints={}):
    """ 
     Checks that :
       - Required parameters are present with the correct type
       - Optional parameters are present, and if not defines them with the default
       - No parameter is Deprecated
       - Constraints are satisfied

     Inputs : 
        - Required   : a dict which associates 
              name -> (documentation string, type_or_type_list) 
        - Optional   : a dict which associates 
              name -> (documentation string,
                       defaultValue : a element of the required type or a function : DIC -> defaultValue,
                       type_or_type_list) 
              where type_or_type_list is a type or a list of types
        - Deprecated : a list of forbidden names for parameters.
        - Constraints : a dict of functions taking the dict and returning bool.

     Action : 
       This method : 
         - checks that all keys of Required are keys of Parameters withthe right type.
         - for all keys of Optional : 
            - if the key is present in DIC, check its type.
            - otherwise, add this key in DIC, with the defaultValue.
         - checks that no parameters has a name in Deprecated.
         - check that every constraints is true when evaluated for DIC. 
     
     Returns : a string containing its report.

       If some type is wrong, or if a required parameters is missing, or a constraints is violated,  
       this function throws an exception.
    """
    report = ""
    
    # First check that parameters are not deprecated
    for p in Deprecated:
        if p in DIC : raise ValueError, "Parameter %s is deprecated."%p

    # Check the Required dict
    for p,v in Required.items() :
        if len(v) not in [2] : raise SyntaxError,"Required : Syntax error for parameter %s"%p

    def check_type(p,val,ty):
       """ check that val is of type ty"""
       if type(ty) not in [ListType,TypeType]:
          raise SyntaxError,"Syntax error for Required parameter %s : second element is not a type or a list of types."%p
       ty_list = ty if type(ty)==ListType else [ty]
       if type(val) not in ty_list : 
            raise TypeError, "%s is not of the correct type. Should be %s, not %s"%(p,ty, type(val))
           
    # Check all required items
    missing={}
    for p,v in Required.items() :    
        if p not in DIC :
            print "missing", p,DIC.keys()
            missing[p] =v
        else :
            doc, ty = v
            check_type(p,DIC[p],ty)

    # If some required are missing, this is a fatal error : list them and exit
    if missing :
        report += "\nI can not find the following required items"
        for p,v in missing.items() :
                report += "\n *- %s of %s \n  Description : %s"%(p,v[1],v[0])
        raise ValueError, report

    # check for the Optional parameters
    report_add, toadd ="\nOptional parameters defined to their default values: ", True
    for p,v in Optional.items() :
        if len(v) == 2 : 
           doc, default= v; ty = type(default)
        elif len(v) ==3 : 
           doc, default, ty = v
        else :
           raise SyntaxError,"Optional : Syntax error for parameter %s"%p
        
        if p in DIC:
           check_type(p,DIC[p],ty)
        else : 
            DIC[p] = default
            if toadd : report, toadd = report + report_add, False
            report += "\n *%s -> %s"%(p,default)
            if doc.strip() : report += "\n   Description :\n     %s\n"%doc
            
    # Check the Constraints
    ok = True
    for name,F in Constraints.items() : 
        if not F(DIC) : 
            ok =False
            report += "\n Constraint '%s' is not satisfied"%(name)
    if not ok : raise ValueError, report
               
    report += "\n------------------------------------------------\n"
    return report 

