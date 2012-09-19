.. index:: dependencies

.. _dependencies:

Dependencies
---------------------


TRIQS is built upon several python and C++ libraries, which, if not present already in your system, can be freely downloaded and installed.
All the libraries and tools used by TRIQS are listed in the table : 

==================    ================  ================================================================================
Libraries/tools       Version           Comment
==================    ================  ================================================================================
mpi                   e.g., openmpi     Parallelism
                                    
                                        Since standard linux distributions (and macports on OS X)
                                        now provides openmpi, even on laptops, we avoid the unnecessary complication
                                        of maintaining a non parallel version of TRIQS
fftw                  >= 3.2            Fourier transform
boost                 >= 1.49           C++ librairies
hdf5                  >= 1.8.0          File storage system. Important: the *serial* version must be installed
python*               >= 2.6.5
scipy*                                  python mathematical library
numpy*                                  python scientific library
h5py*                                   python interface to hdf5 library
sphinx*               >= 1.0.1          python documentation tools
pyparsing*                              Tool for sphinx
matplotlib*           >= 0.99           python 2D plotting library
==================    ================  ================================================================================

 \* designates the libraries included in the Enthought python distribution.

The compilation of TRIQS requires cmake as well as C++ and F90 compilers (the later is used for the Wien2TRIQS interface only).

Tested compilers include : 

* C++

  * g++ 4.4, 4.5 [Linux], 4.2 [Os X]
  * clang++ [ Os X and Linux]
  * icc 11, icc 12 [ Linux]

* F90

  * ifort 
  * gfortran

Quick install of python dependencies using Enthought [Recommended on all platforms]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A simple way to install all the python dependencies is to use the `Enthought <http://www.enthought.com/>`_ scientific python distribution,
which is available for most platforms [tested on version 7.2].

Enthought's distribution allows you to have a recent python, ipython with all the necessary libraries
(including the ipython notebook and its dependencies) without upgrading your machine's system (in particular python).
It is very convenient on clusters, e.g. 

How to proceed : 

* Take the `EPD Free version <http://www.enthought.com/products/epd_free.php>`_  (the free and light version of Enthought).
* Install it ::

   bash epd_free.sh -b -p EPD_Free_Install_Dir

* Add the following packages (the last two are only to compile the documentation) ::

   EPD_Free_Install_Dir/bin/easy_install h5py
   EPD_Free_Install_Dir/bin/easy_install Sphinx
   EPD_Free_Install_Dir/bin/easy_install pyparsing

* When compiling TRIQS, you will simply pass the option ::

    cmake ..... -DPYTHON_INTERPRETER=EPD_Free_Install_Dir/bin/python

 and that is all.

.. warning ::
 
 A priori, you could also use the regular, full Enthought distribution (it is free for academics), but we do not recommend it
 at present. Indeed, in this distribution the HDF5 has been compiled without C++ support, which TRIQS requires.
 So in most machines it will work, with h5py compiled with the Enthought HDF5 version, while the C++ code will silently include
 your system's HDF5 headers, if they exist ...


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



