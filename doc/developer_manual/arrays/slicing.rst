.. highlight:: c

.. _Slicing:

Element access & partial views 
==================================

Element access
--------------------

* The elements of the arrays can be naturally accessed using the ( ) operator::

   array<T,2> A(10,10);               // defines an array
   A(1,2) =10                         // element access.
  
Partial view, slices
------------------------

Various kind of partial views and slices can be made on arrays and matrices.

* A `partial view` is defined as a view of a restricted portion of the array while 
  a `slice` is strictly speaking a partial view of a lower dimension of the original array, 
  e.g. a column of a matrix.

* Partial views uses the ( ) operator, as the evaluation of the array::

   array<T,2> A(10,10);                    // defines an array

   A ( 1, range(0,2) )                     // 1d slice
   A ( 1, range() )                        // 1d slice taking all the second dim
   
   A( range(0,10,2), range(0,10,2))        // a 2d slice viewing every each elements with even coordinates.
   
   array_view<T,1>  SL( A(0,range(0,3)));  // naming the view
   auto SL = A(0,range(0,3));              // C++0x with auto [check this]

* **Return type** : 
  
  * Partial views of array or array_view return an array_view.
  * Partial views of vector or vector_view return an vector_view.
  * 2d partial views of matrix or matrix_view return matrix_view.
  * BUT : (1d) slices of matrix or matrix_view return vector_view.
  * 0d slices of anything are converted to the `value_type` of the array.

The `range` type
^^^^^^^^^^^^^^^^^^^^^

  `range` mimics the python `range`. It can be constructed with : 
  
  * no argument : it then takes the whole set of indices in the dimension (like `:` in python) ::    
    
     A(range(), 0) // take the first column of A

  * two arguments to specify a range ::

     A(range (0,3), 0)  // means  A(0,0), A(1,0), A(2,0)

    .. warning::
       the second element is excluded : range(0,3) is 0,1,2, like in Python. This is different from Blitz++, e.g.

  * three arguments  : a range with a step ::

      A(range(0,4,2), 0) // means A(0,0), A(2,0)  
    

The `ellipsis` type
^^^^^^^^^^^^^^^^^^^^^^

* Ellipsis can be provided in place of `range`, as in python. The type `ellipsis` is similar to range
  except that it is implicitely repeated to as much as necessary.

* Example ::

   array<long,4> B(2,3,4,5) ;
   
   B(0, ellipsis(), 3 )   ; // same as  B(0, range(),range(), 3 ) 
   B(0, ellipsis(),2, 3 ) ; // same as  B(0, range(), 2, 3 ) 
   B( ellipsis(),2, 3 )   ; // same as  B( range(),range(),  2, 3 ) 

* NB : there can be at most one ellipsis per expression (otherwise it would be meaning less).

* Example of usage : 
 
  Ellipsis are useful to write generic algorithms. For example, imagine that you want to sum 
  arrays on their first index ::
        
    template<typename ArrayType>
    array_view<typename ArrayType::value_type, ArrayType::rank-1> sum0 (ArrayType const & A) {
     array<typename ArrayType::value_type, ArrayType::rank-1> res = A(0,ellipsis());
     for (size_t u =1; u< A.shape()[0]; ++u) res += A(u,ellipsis());
     return res;
    }

    ///.....

    array<double,2> A(5,2); 
    array<double,3> B(5,2,3); 
    std::cout<< sum0(A) << sum0(B) <<std::endl;


.. highlight:: python

Bound checking
----------------

To be fast, by default, no check are made on the indices while accessing elements or slicing.
However, checks can be activated in two ways : 

* Adding the `Tag::BoundCheck` option (in the options<...> template of the array, matrix, vector)

* Defining the debug macro TRIQS_ARRAYS_ENFORCE_BOUNDCHECK, which switches the default option from `Tag::NoBoundCheck` to `Tag::BoundCheck`
  for all arrays, matrices and vectors.

In both cases, if the indices are not within the domain of defintion, an exception triqs::arrays::key_error 
will be thrown. It's .what() returns the file and line where the exception occurs, with the stack of all in C++, 
e.g. ::

  BOX>./bound_check_nopy 

  catching key error in B 
  triqs::arrays key_error at ..../triqs_source/triqs/arrays/test/C++/./src/IndexMaps/cuboid/./cuboid_domain.hpp : 104

  Trace is :

  ./bound_check_nopy void triqs::arrays::IndexMaps::cuboid_domain<2>::assert_key_in_domain<boost::tuples::tuple<int, int> > \
  (boost::tuples::tuple<int, int> const&) const 0x181 [0x403e11]
  ./bound_check_nopy main 0x916 [0x403016]
  /lib/libc.so.6 __libc_start_main 0xfd [0x7f389e6abc4d]
  ./bound_check_nopy [0x4023c9]

  key out of domain 
  key [1] = 3 is not within [0,3[

  Further information on the line of the error in the stack can be retrieved with the addr2line utility `[linux only]`
  which retrieve the source file and line from the address of the function (if you compile with -g)  ::

  BOX >addr2line -f -e ./bound_check_nopy 0x403016
  main 
  /scratch/parcolle/triqs_source/triqs/arrays/test/C++/bound_check.cpp:35


