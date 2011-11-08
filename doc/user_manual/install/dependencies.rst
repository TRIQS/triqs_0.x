.. index:: dependencies

.. _dependencies:

Dependencies
---------------------


TRIQS is built upon several python and C++ libraries, which, if not present already in your system, can be freely downloaded and installed.
All the libraries and tools used by TRIQS are listed in the table : 

==================    ============  ================================================================================
Libraries/tools       Version       Comment
==================    ============  ================================================================================
mpi                   openmpi e.g.  Parallelism.
                                    
                                    Since standard linux distributions (and macports on OS X) 
                                    now provides openmpi, even on laptops, we avoid the unnecessary complication
                                    of maintaining a non parallel version of TRIQS.
fftw                  >=3.2         Fourier transform  
boost                 >= 1.46       C++ librairies.
hdf5*                 >= 1.8.x      File storage system
python*               2.6 or 2.7 
scipy*                              python mathematical library         
numpy*                              python scientific library
h5py*                               python interface to hdf5 library
sphinx*               >1.0          Python documentation tools
pyparsing*                          tool for sphinx
matplotlib*           >=0.99        python 2D plotting library
==================    ============  ================================================================================

 '*' designates the libraries included in the Enthought python distribution.

The compilation of TRIQS requires cmake as well as C++ and F90 compilers (the later is used for the Wien2TRIQS interface only).

Tested compilers include : 

* C++

  * g++ 4.4, 4.5 [Linux], 4.2 [Os X]
  * clang++ [ Os X and Linux]
  * icc 11, icc 12 [ Linux]

* F90

  * ifort 
  * gfortran

Quick install of python dependencies using Enthought [all platforms]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A simple way to install all the python dependencies is to use the `Enthought <http://www.enthought.com/>`_ scientific python distribution.

 * It contains recent versions of all the librairies marked with a star in the previous list, i.e. all the python tools and hdf5.  
 * It is free for academic use.
 * It is available for most platforms.

Regular install of dependencies [Ubuntu]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* All the packages and versions have been tested on Ubuntu 10.04 LTS Lucid and higher (with lucid-backports for sphinx version) (except for boost, see above).

Regular install of dependencies [Os X]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* using macports : install the following packages (p26 or py27) ::
      
       sudo port install cmake clang openmpi fftw-3 hdf5-18 py26-matplotlib py26-numpy py26-scipy p26-h5py ... 


Boost [all plaforms]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The boost library is often upgraded, and it is rare to have the correct version in your distribution.  TRIQS installation process offers two choices : 

  * Recommended choice: As explained in the :ref:`page above <installation>`, you can download simply the latest *sources* and TRIQS will do all the job for you by compiling the pieces of boost that are needed in a little boost_for_triqs library.

  * OR you can include and link with an installed boost if the version if high enough as discussed in :ref:`install_options`.



