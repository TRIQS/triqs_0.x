Rational
===========================

.. highlight:: c


This library provides a multi-dimensionnal array library
for numerical computations with the following characteristics/goals : 

* **Simplicity of use**.

  Arrays must be as simple to use as in python (numpy) or fortran.
  This library is designed to be used by physicists, not only by professionnal programmers, 
  We do *a lot* of array manipulations, and we want to maintain *readable* codes.

* **Genericity, abstraction and performance** : 
 
  We want to have simple, readeable code, with the same (or better) speed than manually written low level code.
  Most optimisations should be delegated to the library and the compiler.

* **Complete interoperability with python numpy arrays**.
 
  This library is used a lot with mixed C++/python codes.
  It handles quick conversion between the C++ and python world, e.g. :

   * work on a view of a numpy, 
   * create a array in C++, and return it as a numpy.
   * mix the various kind of arrays transparently in C++ expressions and in cython code.

* **HDF5** : simple interface to hdf5 library for an easy storing/retrieving into/from HDF5 files.

* **Simple interface to (some) blas/lapack** for matrices, vectors.

* **MPI** : compatibility with boost::mpi interface.




