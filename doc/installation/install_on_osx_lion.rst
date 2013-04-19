.. index:: install_on_osx_lion

.. _install_on_osx_lion:
.. highlight:: bash

Installation on Mac OS X [Mountain Lion]
==============================================

Previous versions of Mac OS X are not supported.

On Mountain Lion, clang (llvm) is the default C++ compiler, 
instead of the obsolete gcc 4.2 or previous version.

NB: You can install triqs on previous OS X, but only if you install clang (via Apple) or gcc 4.7 (via macports).

Installation of the dependencies
________________________________

The only supported solution for the mac is homebrew.

#. Install homebrew 
#. Install XCode (directly from the Mac store). In Preferences/Downloads, install "Command Line tools".

#. Install several packages which are needed :: 
         
     brew install cmake
     brew install gfortran
     brew install  --enable-cxx hdf5 
     brew install gsl
     brew install fftw
     brew install open-mpi
     brew install zmq
     brew install python
     brew install doxygen

#. Now install a virtualenv [EXPL. NEEDED] and install python packages ::
    
    pip install h5py
    pip install numpy
    pip install scipy
    pip install git+https://github.com/matplotlib/matplotlib.git#egg=matplotlib-dev
    pip install tornado
    pip install pyzmq
    pip install ipython

#. If you wish to compile the documentation locally, install sphinx, its dependencies and mathjax:: 
  
     pip install sphinx
     easy_install pyparsing==1.5.7
     git clone git://github.com/mathjax/MathJax.git MathJax

   NB : you need pyparsing <1.5.7 since apparently v.2.0 works only for python 3.

#. Download the latest `sources of boost <http://www.boost.org/users/download/>`_  and untar them into a given directory ``BOOST_SRC``


TRIQS installation
__________________

#. Download the TRIQS sources::

      git clone git@github.com:TRIQS/TRIQS.git TRIQS_src

#. Generate a Makefile using cmake::

      cmake TRIQS_src -DBOOST_SOURCE_DIR=BOOST_SRC 

#. Compile TRIQS, its tests and install it into INSTALL_DIR (default) (N is the number of core of your mac)::

      make -jN && make test && make install 

#. If you use Wien2TRIQS, please complete the installation as described :ref:`here <wien2k_inst>`.

