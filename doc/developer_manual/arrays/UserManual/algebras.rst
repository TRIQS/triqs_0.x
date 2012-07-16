
.. highlight:: c

array and matrix/vector algebra 
=======================================================

Arrays and matrices can be combined in formal algebraic expressions, which models the :ref:`HasImmutableArrayInterface` concept.
This algebraic expressions can therefore be used as RHS of assignment (SEE) or in array/matrix contructors.
For example ::
  
   array<long,2> A (2,2), B(2,2),C;
   C= A + 2*B;
   array<long,2> D = A+ 2*B;
   array<double,2> F = 0.5 * A;  // Type promotion is automatic

Expression templates : the basic idea
---------------------------------------------

The programming technique used here is called `expression templates`. 

The principle in a nutshell is that the result of A+B is **NOT** an array, but a complex type of the form::

  auto Expr = A+B; // <--- Expr is a type of the form  addition_of< A_type, B_type>
  auto Expr = A + 2*B + C; // <--- Expr is a type of the form 
  // addition_of<A_type, addition_of< multiplication_of<int,B_type>, C_type> >

The resulting type stores the syntax tree in the type (using the naturally recursive structure of templates)
and it stores the numbers and a view on the arrays involved in the expression.

When an array in assigned (or constructed from) such expression, it fills itself
by evaluating the expression.

This technique allows the elimination of temporaries, so that the clear and readable code::

   Z= A + 2*B + C/2;

is in fact rewritten by the compiler into ::
 
   for (i,j,...) indices of A, B  : 
     C(i,j) = A(i,j) + 2* B(i,j) + C(i,j)/2

instead of making a chain of temporaries (C/2, 2*B, 2*B + C/2...) that "ordinary" object-oriented programming would produce.
As a result, the produced code is as fast as if you were writing the loop yourself,
but with several advantages : 

* It is more **compact** and **readable** : you don't have to write the loop, and the indices range are computed automatically.
* It is much better for **optimization** : 
  
  * What you want is to tell the compiler/library to compute this expression, not *how* to do it optimally on a given machine.
  * For example, since the memory layout is decided at compile time, the library can traverse the data
    in an optimal way, allowing machine-dependent optimization.
  * The library can perform easy optimisations, e.g. for vector it will use blas if possible.
  
Arrays vs matrices
----------------------

Because their multiplication is not the same, arrays and matrices don't form the same algebra.
Mixing them in expression  would therefore be meaningless and it is therefore not allowed ::

   array<long,2> A;
   matrix<long,2> M;
   M + A; // --> ERROR

However, you can always make a matrix_view from a array of rank 2 ::
  
   A + make_matrix_view(M); //--> OK.

.. note::

   Making view is very cheap, it only copies the index systems. 

Expressions are lazy
---------------------------

This means that constructing an expression is separated from evaluating it ::

   auto e =  A + 2*B;             // expression, purely formal, no computation is done
   cout<< e <<endl ;              // prints the expression
   cout<< e(1,2) <<endl ;         // evaluates just at a point
   cout<< e.domain() <<endl ;     // just computes its domain
   array<long,2> D(e);            // now really makes the computation and store the result in D.
   D = 2*A +B;                    // reassign D to the evaluation of the expression.

*N.B.* The expression type is complicated (the expression in stored in the C++ type), so we used here 
the C++11 `auto` to make simple things simple.

FAQ
----------

Where can expressions be used ?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Expressions models :ref:`HasImmutableArrayInterface` concept, so they can used *anywhere* 
an object modeling this concept is accepted, e.g. : 

* array, matrix contruction
* operator =, +=, -=, ...

They behave like an immutable array : they have a domain, they can be evaluated.
When `C` is assigned to the expression in the previous example, 
the compiler just needs to compute the new domain for `C`, resize it and fill it by evaluation of the expression.

What is the cost of this technique ?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Thanks to the Boost Proto library, this can be done :

* with an acceptable increase in compilation time (try it !).
* in about *2 pages of readable code* !

