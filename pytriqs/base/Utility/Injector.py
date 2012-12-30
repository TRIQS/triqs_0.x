
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

import sys,inspect

def make_injector(ClassIntoCodeWillBeInjected) : 
    """
    This is a code injector !

    Usage : 
      If you want to inject code into class A (typically a C++ wrapped class) ::
      
       from pytriqs.base.Utility.Injector import make_injector
       class __inject (make_injector(A) ,Class1, Class 2, A) : # A MUST be the last one
           def method1() :  # etc, as normal class
           def __init__(): 
           # as normal
           self._init_before_injection__(args) # to call the original constructor of A before injection
   
      * After those statements, the methods of Class1, Class2, and __inject will be added to A, **without derivation**.
      * So a C++ code of type ::
         
           return boost::python::object ( an_object_of_the_wrapped_C++_type )
        
        will return an A, with all the python extension.

    .. warning::
       
         Code of Class 1, Class 2 MUST HAVE NO __init__ !


    * doc_processor : reprocess the documentation ...
    """
    class injector(object):
        debug, debug2 = False, False
        class __metaclass__(ClassIntoCodeWillBeInjected.__class__):
            def __init__(self, name, bases, dct):
                if self.debug2 : 
                    print >> sys.stderr, "bases = ", bases
                    print >> sys.stderr, "dct  = ", dct
                target = bases[-1]
                if self.debug2 : print >> sys.stderr, "target", target
                if type(target) not in (self, type):
                    target._init_before_injection__ = target.__init__ # save the original constructor 
                    lst_to_inject = [ b.__dict__ for b in bases[1:-1] ] + [dct]
                    for D in lst_to_inject: 
                        for k,v in D.items():
                            if D!=dct and k=="__init__": 
                                raise KeyError, "no constructor allowed in base classes "
                            if self.debug : print >> sys.stderr, "injecting %s . %s = %s "%(target,k,v)
                            setattr(target,k,v)
                return  type.__init__(self, name, bases, dct)
    return injector


