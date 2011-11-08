.. index:: installation steps

.. _installation:

Installation 
============

Here is a *quick installation* procedure assuming that  
every :ref:`necessary librairies and tools <dependencies>` are present on your machine (see below). 
More specific cases and options are described in the following paragraphs.

#. Get the code::

     git clone git://github.com/TRIQS/TRIQS.git path_to_triqs_sources

#. Download the `sources of boost <http://ipht.cea.fr/triqs/download/boost_1_47_0.tar.bz2>`_ and untar them into ``path_to_boost_sources``.
   **Don't compile boost**.

#. Create a build directory (different from the TRIQS source directory nor a subdirectory thereof) ::
    
    mkdir triqs_build && cd triqs_build 
 
#. Configure ::

    cmake path_to_triqs_sources -DBOOST_SOURCE_DIR=path_to_boost_sources \
          -DCMAKE_INSTALL_PREFIX=path_to_install_directory

   [NB : the install directory is optional, the default value is `triqs_build/INSTALL_DIR`]
 
   Successful run of cmake will finish with ::

     -- Build files have been written to: *path_to_triqs_build*

#. Build TRIQS, test it and install it with (N being the number of cores of your machine) ::
  
    make -jN   && make test && make -jN install 
 
   [NB : the tests SingleSite and CDMFT-v2  take longer than the other ones (30s to 1 min) ].

#. You're done! Now just launch ``pytriqs`` (all the files have been installed under ``path_to_install_directory``).
   
.. toctree::
   :maxdepth: 1

   downloading
   dependencies
   install_options
   ipython 
   wien2k_inst
   where_are_things
   binary

















