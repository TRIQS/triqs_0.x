
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

import sys,re,string
from math import *
import numpy
        
def call_factory_from_dict (cl,dic) :
    """Given a class cl and a dict dic, it calls cl.__factory_from_dict__(dic)"""
    return cl.__factory_from_dict__(dic)

def sum_list(L) :
    """ Can sum any list"""
    return reduce(lambda x, y: x+y, L) if len(L)>0 else []

def conjugate(x) :
    try:
      return x.conjugate()
    except:
      print "The object you're using conjugate on has no attribute conjugate!"

def real(x) : 
     return x.real if hasattr(x,'real') else x

def imag(x) : 
     return x.imag if hasattr(x,'imag') else x

def sign(x):
    if x>0.0 : return 1
    if x<0.0 : return -1
    return 0

def cartesian_product(*L):
    """The cartesian product of lists"""
    return [li for li in x_cartesian_product(*L)]

def x_cartesian_product(*L):
    """ A generator for the cartesian product of lists"""
    if len(L)==1 :
        for x in L[0] :
            yield tuple(x)
    if len(L)==2 :
        for x in L[0]:
            for y in L[1]:
                yield (x,y)
    else:
        for x in L[0]:
            for xx in apply(x_cartesian_product,L[1:]):
                yield tuple([x] + list(xx))


def mustbe(test,err,mess):
    if not(test) :
        raise err,mess
def mustnotbe(test,err,mess):mustbe(not(test),s,mess)

def is_array(x): return(type(x)== type(numpy.array([1.0])))
def is_dict(x): return(type(x)== type({}))
def is_string(x): return(type(x)== type(""))
def is_type(x,l): return(type(x) in map(lambda x : type(x),l))
#def is_list(x): return is_type(x,[ [] ])
#def is_list_or_tuple(x): return is_type(x,[ [],() ])
def is_number(x): return(type(x) in [type(1), type(1j), type(1.0)])


def nint_strict(x,precision=1.e-9):
    """ Round x to the closest integer and asserts that its distance to this integer is less than precision.
        precision must satisfy :  precision >0 and precision <0.5
    """
    assert precision >0 and precision <0.5, "nint_strict : precision makes no sense !"
    if is_array(x) :
        i = numpy.floor(x+0.5)
        #print i, abs(i-x), max(abs(i-x))
        #assert abs(i-x).max()<precision,  repr(i) + "\n "+repr(x) + "\n The Float is not close enough to the integer "
        assert abs(i-x).max() <precision,  repr(i) + "\n "+repr(x) + "\n The Float is not close enough to the integer "
        return i.astype(numpy.int)
    else :
        i = floor(x+0.5)
        assert abs(i-x)<precision,  repr(i) + " "+repr(x) + " The Float is not close enough to the integer "
        return (int(i))
    
class PrettyPrint:
    def __init__(self, exclude=()) :
        self._excl = exclude
        
    def __repr__(self):
        """ Print all the elements of the class which are list, string, number, bool, Numarrays"""
        _type_to_print = [type([]),type(1),type(1.0),type(True),type(numpy.array(1)),type({}),type('')]
        r=self.__class__.__name__ + ' class : \n\n   '
        for nom,val in self.__dict__.items() + [ x for x in self.__class__.__dict__.items() if x not in self.__dict__.items()] :
            if (type(val) in _type_to_print ) and  nom[0]!='_' and nom not in self._excl :
                r= r+ nom + " : " + repr(val) + "\n   "
        for nom in self._excl :
            r= r+ nom + " : Not shown" + "\n   "

        return(r)

def find_instance(CLASS,DIC):
    """ Return all instance of the CLASS in the dictionnay DIC"""
    return([ y for x,y in DIC.items() if isinstance(y,CLASS)])


def s_to_number(x, return_type):
    """
    Transform a string into a complex or a float (using a regexp).
    x : string containing the number
    return_type : 'float' or 'complex'
    Returns : if one  number has been found, this number, otherwise a list of the number found.
    """
    def assemble(t) : 
       return float(t[0]) * (  10** float(t[1])  if t[1] else 1)

    if type(x) !=type('') : raise TypeError, "x must be a string"
    if return_type in ['float','complex'] :
      reg ={'float' : r'([-+]?\d+(?:\.\d*)?|\d*\.\d+)(?:[eEdD]([-+]?\d+))?'}
      reg['complex'] =  r"\(\s*%s\s*,\s*%s\s*\)"%(reg['float'],reg['float']) 
      res1 = re.compile(reg[return_type]).findall(x)
      #print res1
      res = []
      if return_type=='float':
        for t in res1 :
          res.append( assemble(t[0:2]))
      else:
        for t in res1 :
          res.append(complex(assemble(t[0:2]),assemble(t[2:4])))
    else :
      raise TypeError, "return_type unknown"
    if len(res)==1 : return res[0]
    else : return res

if __name__=='__main__':
    print s_to_number('-3.0596071669559', 'float')
    print s_to_number('9.1e-05', 'float')
    print s_to_number('9.1e05', 'float')
    print s_to_number('(-3.0596071669557,0)','float')
    print s_to_number('(-3.0596071669557,2.e-3)','complex')

