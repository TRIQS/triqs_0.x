
Reference manual
##############################

The :class:`HDF_Archive` class offers a convenient interface between the python objects and the HDF5 files, 
similar to a dictionary (or a shelve).

The module contains two classes : 

* :class:`HDF_Archive_group` : operates on a subtree of the HDF5 file
* :class:`HDF_Archive` : an :class:`HDF_Archive_group` at the root path, with a simple constructor.
* :class:`HDF_Archive_Inert` : 

Typically, one constructs an :class:`HDF_Archive` explicitely, the :class:`HDF_Archive_group` is created during operations, e.g.::

  h = HDF_Archive( "myfile.h5", 'r') 
  g = h['subgroup1'] # g is a HDF_Archive_group.
  
Apart from the root path and the constructor, the classes are the same (in fact :class:`HDF_Archive` `is a` :class:`HDF_Archive_group`).
Let us first document :class:`HDF_Archive`.

.. warning::  

   :class:`HDF_Archive` and :class:`HDF_Archive_group` do **NOT** handle parallelism.
   Check however the :class:`HDF_Archive_Inert` below.   
 
HDF_Archive
=======================================

.. autoclass::  pytriqs.base.archive.HDF_Archive



HDF_Archive_group
=======================================

.. class::  HDF_Archive_group

   There is no explicit constructor for the user of the class.

   The :class:`HDF_Archive_group` support most of the operations supported by dictionaries. In the following, *H* is a :class:`HDF_Archive_group`. 

   .. describe:: len(H)

      Return the number of items in the :class:`HDF_Archive_group` *H*.

   .. describe:: H[key]

      Return the item of *H* with key *key*, retrieved from the file.  Raises a :exc:`KeyError` if *key*
      is not in the :class:`HDF_Archive_group`.

   .. method:: get_raw (key)

      Returns the subgroup key, **without any reconstruction**, ignoring the HDF5_data_scheme.

   .. describe:: H[key] = value

      Set ``H[key]`` to *value*.

   .. describe:: del H[key]

      Remove ``H[key]`` from *H*.  Raises a :exc:`KeyError` if *key* is not in the
      :class:`HDF_Archive_group`.

   .. describe:: key in H

      Return ``True`` if *H* has a key *key*, else ``False``.

   .. describe:: key not in H

      Equivalent to ``not key in H``.

   .. describe:: iter(H)

      Return an iterator over the keys of the dictionary.  This is a shortcut
      for :meth:`iterkeys`.

   .. method:: items()

      Generator returning couples (key, values) in the group.

   .. warning::

      Note that in all these iterators, the objects will only be retrieved from the file and loaded into memory 
      one by one. 


   .. method:: keys()

      Generator returning the keys of the group.

   .. method:: update(d)

      Add into the archive the content of any mapping *d* : keys->values, with 
      hfd-compliant values.

   .. method:: values()

      Generator returning the values in the group
      
   .. method:: create_group(K)

      Creates a new subgroup named `K` to the root path of the group. 
      Raises exception if the subgroup already exists.
  
   .. method:: is_group(K)
   
      Return True iif K is a subgroup.

   .. method::  is_data(K)
   
      Return True iif K is a leaf.

   .. method:: read_attr(AttributeName) 

      Return the attribute `AttributeName` of the root path of the group.
      If there is no attribute, return None.

   .. method:: root_path()

      Return the root path of the group

   .. method:: apply_on_leaves (f) 

       For each named leaf (name,value) of the tree, it calls `f(name,value)`.

       f should return : 
      
           - `None`                    : no action is taken
           - an `empty tuple` ()       : the leaf is removed from the tree
           - an hdf-compliant value    : the leaf is replaced by the value
 
HDF_Archive_Inert
=======================================


.. class:: HDF_Archive_Inert
 
   :class:`HDFArchive` and :class:`HDFArchive_group` do **NOT** handle parallelism.
   In general, one wish to write/read only on master node, which is a good practice
   cluster : reading from all nodes may lead to communication problems.

   To simplify the writing of code, the simple HDF_Archive_Inert class may be useful.
   It is basically inert but does not fail.

   .. describe:: H[key]

      Return H and never raise exception so e.g. H['a']['b']  never raise exception...

   .. describe:: H[key] = value
       
      Does nothing

   Usage in a mpi code, e.g. ::

     R = HDF_Archive("Results.h5",'w') if MPI.IS_MASTER_NODE() else HDF_Archive_Inert()
     a= MPI.bcast(R['a'])       # properly broadcast the R['a'] from the master to the nodes.
     R['b'] = X                 # sets R['b'] in the file on the master only, does nothing on the nodes.


.. highlight:: python

.. _HDF_Protocol:

Hdf-compliant objects
=======================================

By definition, hdf-compliant objects are those which can be stored/retrieved in an :class:`HDF_Archive`.
In order to be hdf-compliant, a class must : 

* have a `HDF5_data_scheme` tag properly registered.
* implement one of the two protocols described below.

.. _HDF_Protocol_data_scheme:

HDF5 data scheme
-----------------------

To each hdf-compliant object, we associate
a `data scheme` which describes how the data is stored in the hdf5 tree, 
i.e. the tree structure with the name of the nodes and their contents.
This `data scheme` is added in the attribute `HDF5_data_scheme` 
at the node corresponding to the object in the file.

For a given class `Cls`, the `HDF5_data_scheme` is `Cls._hdf5_data_scheme_` if it exists or the name of the class `Cls.__name__`.
The `HDF5_data_scheme` of a class must be registered in order for :class:`HDF_Archive` to properly reconstruct the object when rereading.
The class is registered using the module `HDF_Archive_Schemes` ::

 class myclass :
   pass #....
   
 from pytriqs.base.archive.hdf_archive_schemes import register_class
 register_class (myclass)

The function is 

.. py:function:: register_class (cls [, doc = None] ) 

     :param cls: the class to be registered.
     :param doc: a doc directory

     Register the class for :class:`HDF_Archive` use.
      
     The name of `data scheme` will be `myclass._hdf5_data_scheme_` if it is defined, or the name of the class otherwise.

.. _HDF_Protocol_details:

How to become hdf-compliant ?
-----------------------------------

There are two ways in which a class can become hdf-compliant: 

.. toctree::
  :maxdepth: 1

  protocol1
  protocol2

