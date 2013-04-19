.. index::
  single: Green's functions; full Green's function
  module: gf.local

.. _fullgreen:

The complete Green's function (BlockGf) 
=======================================

As mentioned in the introduction, due to the symmetry, a local Green's function usually
has a block structure.
We refer to this object as the `full` or `complete` Green's function, in contrast to
the blocks it is made of.

Most properties of this object can be remembered by the simple sentence:

`A full Green's function is an ordered dictionary name -> block, or equivalently a list of tuples (name, block).`

The blocks can be any of the matrix-valued Green's functions described :ref:`above<blockgreen>`.
The role of this object is to gather them, and simplify the code writing
by factorizing some simple operations.

A little example
--------------------

To start with an example, imagine that the problem
that we consider could involve 5 d-bands of a solid that, for symmetry reasons,
are separated into 2 eg and 3 t2g bands. We therefore first construct the 2
corresponding block Green's functions (in Matsubara frequencies for example)
and group these blocks into a full Green's function `G` with::

  from pytriqs.gf.local import *
  g1 = GfImFreq(indices = ['eg1','eg2'], beta = 50, n_points = 1000, name = "egBlock")
  g2 = GfImFreq(indices = ['t2g1','t2g2','t2g3'], beta = 50, n_points = 1000, name = "t2gBlock")
  G = BlockGf(name_list = ('eg','t2g'), block_list = (g1,g2), make_copies = False)

where:

* `NameList` is the ordered list of the names of the blocks. 
* `BlockList` is the corresponding list of block Green's function. 
* `Copy` lets you specify if the blocks of the full Green's function are **copies** of the blocks given in
* `BlockList` of if they are **views** of these blocks, see :ref:`below<fullgreencopypolicy>`

These names will be used when we try to access a particular block, for example ::

  >>> G
  Green's Function  composed of 2 blocks at inverse temperature Beta = 50.0: 
  GfImFreq eg :  Beta = 50.000; IndicesL = ['eg1', 'eg2'], IndicesR = ['eg1', 'eg2']  
  GfImFreq t2g :  Beta = 50.000; IndicesL = ['t2g1', 't2g2', 't2g3'], IndicesR = ['t2g1', 't2g2', 't2g3'] 
  
  >>> G['eg']
  GfImFreq eg :  Beta = 50.000; IndicesL = ['eg1', 'eg2'], IndicesR = ['eg1', 'eg2'] 


Reference 
----------------

.. autoclass:: pytriqs.gf.local.BlockGf
  :members: copy, copy_from
 

Operations
---------------

The full Green's functions support various simple operations, that are simply done block by block.

.. note::
   All these operations compute the array of data, but also, if present in the object, the high frequency expansion tail automatically.

* compound operators, `+=`, `-=`, `*=`, `\=` : RHS can be a Green's function of the same type or an expression

* arithmetic operations : `+`, `-`, `*`, `/`, e.g. ::
   
   G = G1 + 2*G2

* inversion, e.g. ::
  
   inv = inverse(g)
   g2=inverse(inverse( g ) - sigma) # this is a Dyson equation

Block access
----------------

Blocks can be accessed like in a `dict` :
  
These names will be used when we try to access a particular block, for example ::
    
  G['eg']

The generic way to access a Green's function element :math:`G^a_{i j}` is therefore ::
  
  G[a][i,j]

Iterator
--------------------

One can iterate on the blocks ::

    for name, g in G: 
      do_something()

In the example above ::
  
  >>> for name, g in G:
  ...  print name, g
  
  eg GfImFreq eg :  Beta = 50.000; IndicesL = ['eg1', 'eg2'], IndicesR = ['eg1', 'eg2'] 
  t2g GfImFreq t2g :  Beta = 50.000; IndicesL = ['t2g1', 't2g2', 't2g3'], IndicesR = ['t2g1', 't2g2', 't2g3'] 


As a result ::
        
    BlockGf( name_block_generator= G, copy=False) 
    
generates a new Green's function `G`, viewing the same blocks.
More interestingly ::
      
    BlockGf( name_block_generator= [ (index,g) for (index,g) in G if Test(index), copy=False)]


makes a partial view of some of the blocks selected by the `Test` condition.

.. warning::
   The order in which the blocks appear is guaranteed to be the same as in the constructor.
   This is why the Green's function is similar to an **ordered** dictionary, not a simple dict.

View or copies ?
---------------------
.. _fullgreencopypolicy:

The Green's function is to be thought like a dict, hence accessing the 
block returns references. When constructing the Green's function BlockGf, 
the parameter `Copy` tells whether a copy of the block must be made before 
putting them in the Green function or not.

.. note::
   This is the standard behaviour in python for a list of a dict.

Example: 

* If you define a Green's function with::

   G = BlockGf(name_list = ('eg','t2g'), block_list = (g1,g2), make_copies = False)

  .. note:: 
 
    Copy is optional, False is the default value. We keep it here for clarity.

  The ``Copy = False`` implies that the blocks of ``G`` are *references* ``g1`` and ``g2``.
  So, if you modify ``g1``, say by putting it to zero with ``g1.zero()``, then the
  first block of G will also be put to zero. Similarly, imagine you define two
  Green's functions like this::

   G1 = BlockGf(name_list = ('eg','t2g'), block_list = (g1,g2), make_copies = False)
   G2 = BlockGf(name_list = ('eg','t2g'), block_list = (g1,g2), make_copies = False)

  Here G1 and G2 are exactly the same object, because they both have blocks
  which are views of ``g1`` and ``g2``. 

* Instead, if you write::

   G = BlockGf(name_list = ('eg','t2g'), block_list = (g1,g2), make_copies = True)

  The ``Copy = True`` ensures that the blocks of G are new copies of ``g1``
  and ``g2``. If you then modify ``g1`` it will not have any effect on G.
  Clearly if you define::

   G1 = BlockGf(name_list = ('eg','t2g'), block_list = (g1,g2), make_copies = True)
   G2 = BlockGf(name_list = ('eg','t2g'), block_list = (g1,g2), make_copies = True)

  Here ``G1`` and ``G2`` are different objects, both having made copies
  of ``g1`` and ``g2`` for their blocks.

  An equivalent writing is ::

    G1 = BlockGf(name_list = ('eg','t2g'), block_list = (g1.copy(),g2.copy()))
    G2 = BlockGf(name_list = ('eg','t2g'), block_list = (g1.copy(),g2.copy()))

shelve / pickle 
---------------------

Green's functions are `pickable`, i.e. they support the standard python serialization techniques.

* It can be used with the `shelve <http://docs.python.org/library/shelve.html>`_ and `pickle <http://docs.python.org/library/pickle.html>`_  module::
  
     import shelve
     s = shelve.open('myfile','w')
     s['G'] = G  # G is stored in the file.

* It can be sent/broadcasted/reduced over mpi ::

     from pytriqs.utility import MPI
     mpi.send (G, destination)

.. warning::
   Shelve is not a portable format, it may change from python version to another (and it does).
   For portability, we recommend using the HDF5 interface for storing data on disks.


HDF5
--------

BlockGf are hdf-compatible with the following HDF5 data scheme

The BlockGf(TRIQS_HDF5_data_scheme = "BlockGf") is decomposed in the following objects : 

=========================   ===========================  ===========================================================================
Name                        Type                         Meaning
=========================   ===========================  ===========================================================================
__BlockIndicesList          string                       The python repr of the list of blocks, e.g.  ('up', 'down')
__Name                      string                       Name of the Green's function block
__Note                      string                       Note 
For each block name         type of the block            The Block Green's function 
=========================   ===========================  ===========================================================================

Example::

 /Gtau                    Group
 /Gtau/__BlockIndicesList Dataset {SCALAR}
 /Gtau/__Name             Dataset {SCALAR}
 /Gtau/__Note             Dataset {SCALAR}
 /Gtau/down               Group
 /Gtau/down/Data          Dataset {1, 1, 1000}
 /Gtau/down/ ... 
 ...
 /Gtau/up                 Group
 /Gtau/up/Data            Dataset {1, 1, 1000}
 /Gtau/up/ ...



