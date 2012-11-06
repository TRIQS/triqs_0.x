.. highlight:: c

Automatic assignment 
=======================

It is possible to use a subset of possible expressions as Left Hand Side (LHS) in an assignment statement, e.g. ::

 A(x_)        = some_expression_of_x_
 A[i_]        = some_expression_of_i_
 A(x_)(i_,j_) = some_expression_of_x_i_j
 A[x_](i_,j_) = some_expression_of_x_i_j


* Right Hand Side (RHS) of the = statement can be any expression.
* Left Hand Side (LHS) of the = sign. A must be a `lazy-assignable`, called by [] or (), one or several time.

This writing means that a regular function ::
  
  x_ -> some_expression_of_x_

will be given to A so that it can fill itself by evaluating this function.
A determines the range of value of x on which the function is called.

How this is done in details depends on the object.


Object with HasRandomAccess concept (i.e. std::vector like)
------------------------------------------------------------

The HasRandomAccess concept is a subset of std::vector ::
  
  value_type & operator[](size_t);
  size_t size() const; 
  
If A models HasRandomAccess then ::

 A[i_] = some_expression_of_i

is rewritten by the triqs::clef lib into ::

 for (size_t u=0; u<A.size(); ++u) A[u] = eval(some_expression_of_i, i_=u)

In practice, for such an object T the trait ::

   triqs::clef::has_random_access_concept<T>

must be true.

In practice : 

* A std::vector models this concept. Just include ::

  #include <triqs/clef/adapters/vector.hpp>
  
* For a custom type, defines e.g. ::

   namespace triqs { namespace clef { 
     template<> struct has_random_access_concept< My_Type >: mpl::true_{};  
   }}
 

.. _callable_object:

Object with a set_from_function method
--------------------------------------------

For objects that have a set_from_function method::

 F(x_,y_) = an_expression_of(x_,y_);

is rewritten by the library into ::

 F.set_from_function( triqs::clef::make_function( an_expression_of(x_,y_), x_, y_) );

In practice, the object must implement a (template) method void set_from_function(), which takes a function and uses
it to recompute the data of the object. The library will autodetect whether this function is present and use it,
unless the trait ::

   namespace triqs { namespace clef { 
    template <typename RHS> struct has_set_from_function<The_Type,RHS>; 
   }}
is specialized to false (where The_Type is the type of the object).

Example ::

  struct MyObject { 
   
  // ....


  template<typename Fnt> void set_from_function (Fnt f) { 
     // do something with the function f
   }
  };
  
  //---------------------------------------

  placeholder <1> x_;
  placeholder <2> y_;
  MyObject ob;
  ob(x_,y_) = 2*x_ + 3*y_; 
  // same as
  ob.set_from_function( make_function(  2*x_ + 3*y_, x_,y_));

General case 
-----------------------------

To be written. Specializing the assignment_impl trait.

