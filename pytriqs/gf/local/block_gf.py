
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

from itertools import izip
import operator
from impl_plot import PlotWrapperPartialReduce

def call_factory_from_dict (cl,dic) :
    """Given a class cl and a dict dic, it calls cl.__factory_from_dict__(dic)"""
    return cl.__factory_from_dict__(dic)

class BlockGf(object):
    """
    Generic Green Function by Block. 
    """
    # not implemented as a dictionnary since we want to be sure of the order !).
    def __init__(self, **kwargs) : 
        """
   * There are several possible constructors, which accept only keyword arguments.
            
            * BlockGf(name_list = list of names, block_list = list of blocks, make_copies = False, name = '')
                   
                   * ``name_list`` : list of the name of the blocks : e.g. ["up","down"].
                   * ``block_list`` : list of blocks of Green functions.
                   * ``make_copies`` : If True, it makes a copy of the blocks and build the Green function from these copies.
            
            * BlockGf(name_block_generator, make_copies = False, name = '')
                   
                   * ``name_block_generator`` : a generator of (index, block)
                   * ``make_copies`` : If True, it makes a copy of the blocks and build the Green function from these copies.
                     
          """
        # first extract the optional name argument
        self.name = kwargs.pop('name','G')
        self.note = kwargs.pop('note','')
        self._rename_gf = kwargs.pop('rename_gf',True)

        if 'make_copies' not in kwargs : kwargs['make_copies'] = False
 
        if set(kwargs.keys()) == set(['name_list','block_list','make_copies']) :
            BlockNameList, GFlist = kwargs['name_list'],kwargs['block_list']
        elif set(kwargs.keys()) == set(['name_block_generator','make_copies']) : 
            BlockNameList,GFlist = zip (* kwargs['name_block_generator'])
        else : 
            raise RuntimeError, "BlockGf construction : error in parameters, see the documentation"
        if kwargs['make_copies'] : GFlist =[g.copy() for g in GFlist]

        # First a few checks
        assert GFlist !=[], "Empty list of blocks !"
        for ind in BlockNameList : assert str(ind)[0:2] !='__', "indices should not start with __"
        assert len(set(BlockNameList)) == len(BlockNameList),"Bloc indices of the Green Function are not unique"
        assert len(BlockNameList) == len(GFlist), "Number of indices and of Green Function Blocks differ"

        # All blocks are compatible for binary operation
        # --> correction : All blocks have the same type
        #if not reduce (operator.and_,[ GFlist[0]._is_compatible_for_ops(x) for x in GFlist[1:] ] , True) :
        #    raise RuntimeError, "The blocks are not compatible for binary operations : not the same type, same temperature, etc..."
        if len(set([ type(g) for g in GFlist])) != 1 :
            raise RuntimeError, "BlockGf : All block must have the same type %s"%GFlist

        # init
        self.__indices,self.__GFlist = BlockNameList,GFlist
        try : 
            self.__me_as_dict = dict(self)
        except TypeError: 
            raise TypeError, "indices are not of the correct type"
        self.__BlockIndexNumberTable = dict( (i,n) for n,i in enumerate(self.__indices) ) # a dict : index -> numero of its order
        
        # Add the name to the G
        self.note = ''
        if self._rename_gf:
            for i,g in self : g.name = "%s_%s"%(self.name,i) if self.name else '%s'%(i,)
        del self._rename_gf

    #------------ copy and construction -----------------------------------------------

    def copy(self,*args):
       """Returns a (deep) copy of self (i.e. new independant Blocks initialised with self) """
       return self.__class__ (name_list = self.__indices[:], block_list = [ g.copy(*args) for g in self.__GFlist],make_copies=False)
    
    def view_selected_blocks(self, selected_blocks) :
        """Returns a VIEW of the selected blocks of self, in the same order as self."""
        for b in selected_blocks : assert b in self.__indices,"Selected Blocks must be existing blocks"
        return self.__class__ ( name_block_generator = [(n,g) for n,g in G if n in selected_blocks ],make_copies=False)

    def copy_selected_blocks(self, selected_blocks) :
        """Returns a COPY of the selected blocks of self, in the same order as self."""
        for b in selected_blocks : assert b in self.__indices,"Selected Blocks must be existing blocks"
        return self.__class__ ( name_block_generator = [(n,g) for n,g in G if n in selected_blocks ],make_copies=True)
 
    def copy_from(self, G2):
        """Copy the green function from G2: G2 MUST have the same structure !!""" 
        assert isinstance(G2, BlockGf)
        for (i,g),(i2,g2) in itertools.izip(self,G2) : 
           if  (g.N1,g.N2) != (g2.N1,g2.N2) : 
               raise RuntimeError, "Blocks %s and %s of the Green Function do have the same dimension"%(i1,i2) 
        for (i,g),(i2,g2) in itertools.izip(self,G2) : g.copy_from(g2)

     #--------------  Iterators -------------------------

    def __iter__(self) : 
        return izip(self.__indices,self.__GFlist)

    #---------------------------------------------------------------------------------

    def _first(self) :
        return self.__GFlist[0]
        
    def __len__(self) : 
        return len(self.__GFlist)

    #---------------- Properties -------------------------------------------

    # Deprecated. Left for backward compatibility ?? but  
    # it would not work for Legendre 
    @property
    def mesh(self) : 
        """Deprecated : backward compatibility only"""
        return  self._first().mesh

    @property
    def beta(self) : 
        """ Inverse Temperature"""
        return  self._first().beta
    
    @property
    def indices(self) : 
        """A generator of the block indices"""
        for ind in self.__indices: 
            yield ind
    
    @property 
    def all_indices(self):
       """  An Iterator on BlockGf indices and indices of the blocs of the form : block_index,n1,n2, where n1,n2 are indices of the block"""
       for sig,g in self : 
          val = g.indices
          for x in val :
              for y in val : 
                  yield  (sig,x,y)

    @property 
    def n_blocks(self):
        """ Number of blocks"""
        return len(self.__GFlist)

    #----------------------   IO    -------------------------------------
    
    def __mymakestring(self,x):
        return str(x).replace(' ','').replace('(','').replace(')','').replace("'",'').replace(',','-')
    
    def save(self, filename, accumulate=False) : 
        """ Save the Green function in i omega_n (as 2 columns).
           - accumulate : 
        """
        for i,g in self : 
            g.save( "%s_%s"%(filename, self.__mymakestring(i)), accumulate)
 
    def load(self, filename, no_exception = False):
        """ 
        adjust_temperature : if true, 
        """
        for i,g in self : 
            try : 
                g.load( "%s_%s"%(filename,self.__mymakestring(i)))
            except : 
                if not(no_exception) : raise  
    
    def __reduce__(self):
        return call_factory_from_dict, (self.__class__,self.__reduce_to_dict__())

    def __reduce_to_dict__(self):
        val = {'name' : self.name, "note": self.note, "indices" : repr(self.__indices) }
        val.update( dict(self ) )
        return val 

    @classmethod
    def __factory_from_dict__(cls,value):
        exec "ind=%s"%value['indices']
        res = cls(name_list = ind, block_list = [value[i] for i in ind], make_copies=False, name=value['name'])
        res.note = value['note']
        return res

    #--------------  Pretty print -------------------------

    def __repr__(self) :
        s =  "Green Function %s composed of %d blocks : \n"%(self.name,self.n_blocks)
        #s =  "Green Function %s composed of %d blocks at inverse temperature beta = %s: \n"%(self.name,self.n_blocks,self.beta)
        for i,g in self:
            s += " %s \n"%repr(g)  #"  Bloc %s  : indices %s \n"%(i,self[i].indices)
        if self.note : s += "NB : %s\n"%self.note
        return s

    def __str__ (self) : 
           return self.name if self.name else repr(self)
 
    #--------------  Bracket operator []  -------------------------
    
    def __getitem__(self,key):
        try :
            g = self.__me_as_dict[key]
        except KeyError : 
            raise IndexError, "bloc index '" + repr(key) + "' incorrect. Possible indices are : "+ repr(self.__indices)
        return g

    def __setitem__(self,key,val):
        try :
            g = self.__me_as_dict[key]
        except KeyError : 
            raise IndexError, "bloc index '" + repr(key) + "' incorrect. Possible indices are : "+ repr(self.__indices)
        g <<= val
       
    # -------------- Various operations -------------------------------------

    def __le__(self, other) : 
        raise RuntimeError, " Operator <= not defined "

    def __ilshift__(self, A): 
        """ A can be 2 things :
          * G <<= any_init will init all the BlockGf with the initializer
          * G <<= g2 where g2 is a BlockGf will copy g2 into self
          """
        if isinstance(A, self.__class__):
           for (i,g) in self : g.copy_from(A[i])
        else: 
           for i,g in self:  g <<= A
        return self
   
    def __iadd__(self,arg):
        if isinstance(arg, self.__class__):
            for (i,g),(i2,g2) in izip(self,arg) : g += g2
        elif operator.isSequenceType(arg):
            assert len(arg) == len(self.__GFlist), "list of incorrect length"
            for l,g in izip(arg,self.__GFlist) : g +=l 
        else:
            for i,g in self: g += arg
        return self

    def __add__(self,y):
        c = self.copy()
        c += y
        return c

    def __radd__(self,y): return self.__add__(y)

    def __isub__(self,arg):
        if isinstance(arg, self.__class__) :
           for (i,g),(i2,g2) in izip(self,arg) : g -= g2
        elif operator.isSequenceType(arg):
            assert len(arg) == len(self.__GFlist) , "list of incorrect length"
            for l,g in izip(arg,self.__GFlist) : g -=l 
        else:
            for i,g in self: g -= arg
        return self

    def __sub__(self,y):
        c = self.copy()
        c -= y
        return c

    def __rsub__(self,y):
        c = (-1)*self.copy()
        c += y
        return c

    def __imul__(self,arg):
        if isinstance(arg, BlockGf): 
            for (i,g),(i2,g2) in izip(self,arg) : g *= g2
        elif operator.isSequenceType(arg) :
            assert len(arg) == len(self.__GFlist) , "list of incorrect length"
            for l,g in izip(arg,self.__GFlist) : g*=l 
        else : 
            for i,g in self : g *= arg
        return self

    def __mul__(self,y):
        c = self.copy()
        c *= y
        return c

    def __rmul__(self,x): return self.__mul__(x)

    def __idiv__(self,arg):
        if operator.isSequenceType(arg) :
            assert len(arg) == len(self.__GFlist) , "list of incorrect length"
            for l,g in izip(arg,self.__GFlist) : g /=l 
        else : 
            for i,g in self : g /= arg
        return self

    def __div__(self,y):
        c = self.copy()
        c /= y
        return c


  
   #-----------------------------plot protocol -----------------------------------

    @property
    def real(self) : 
        """Use self.real in a plot to plot only the real part"""
        return PlotWrapperPartialReduce(self, RI='R')

    @property
    def imag(self) :
        """Use self.imag in a plot to plot only the imag part"""
        return PlotWrapperPartialReduce(self, RI='I')

    def _plot_(self, *l, **kw) : 
        """ Implement the plot protocol"""
        # NB : it is important to make a copy of the dictionnary, since the
        # options are consumed by the _plot_ method
        # however, a shallow copy should be enough here, no need for a deep copy
        return sum([ g._plot_(*l, ** dict (** kw)) for sig,g in self],[])

   #--------------------------------------------------------------------------

    def zero(self) : 
        for i,g in self : g.zero()
  
    def density(self) : 
        """Returns the density as a dict of density of the blocks"""
        return dict( (s,g.density()) for s,g in self )

    def total_density(self) : 
        """ Total density of G  """
        return sum([ g.total_density()  for i,g in self ]).real

    def __check_attr(self,ATTR) :
        if not hasattr(self._first(), ATTR ) :
           raise RuntimeError, "The blocks of this Green Function do not possess the %s method"%ATTR

    def invert(self) : 
       """Inverse all the blocks"""
       self.__check_attr("invert")
       for i,g in self : g.invert()

    def delta(self) :
       """Compute delta from G0"""
       self.__check_attr("delta")
       return self.__class__( name_block_generator = [ (n, g.delta()) for n,g in self], make_copies=False)

    def transpose(self):
       """Transpose of the BlockGf"""
       self.__check_attr("transpose")
       return self.__class__( name_block_generator = [ (n, g.transpose()) for n,g in self], make_copies=False)

    def conjugate(self):
       """Conjugate of the BlockGf"""
       self.__check_attr("conjugate")
       return self.__class__( name_block_generator = [ (n, g.conjugate()) for n,g in self], make_copies=False)

#---------------------------------------------------------

from pytriqs.archive.hdf_archive_schemes import register_class
register_class (BlockGf)

