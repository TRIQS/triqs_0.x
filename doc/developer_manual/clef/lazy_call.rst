
.. highlight:: c

.. _callable_object:

Overloading functions and callable objects for lazy expression arguments 
=============================================================================

Given a callable object or a function, it is possible to overload it for lazy expression arguments, returning a lazy expression.
In all  cases:

   *Overloads that are activated iif at least one argument is a lazy expression*.

Functions
------------

The `TRIQS_LAZY_MAKE_FNT_LAZY` macro defines the overload for lazy expressions arguments of a function. Synopsis is ::

 namespace triqs { namespace lazy { 
 TRIQS_LAZY_MAKE_FNT_LAZY (function_arity, function_to_make_lazy);
 }}

For example:: 

 #include <triqs/clef/core.hpp>
 double foo(double x) { return x/2;}
 int foo(int x) { return x/2;}

 namespace triqs { namespace clef { 
  TRIQS_LAZY_MAKE_FNT_LAZY (1, foo) ;
 }}

Note that : 
 
 * This overload **must** be defined in the triqs::clef namespace, since it is found by ADL.
 * The function `foo` can have many overloads.
 * The function foo can be a template, BUT then the template must be disabled for lazy expressions, as in ::

    #include <triqs/clef/core.hpp>
       
    template<typename T> 
    typename boost::disable_if<triqs::clef::is_lazy<T>,T>::type 
    foo (T const & x) { return x+1;}
      
    namespace triqs { namespace clef { 
     TRIQS_LAZY_MAKE_FNT_LAZY (1, foo) ;
    }}



Callable objects
--------------------

Similarly to functions, classes can define an `operator()` for lazy expressions arguments.
This is done with one of the two macros :

* TRIQS_LAZY_ADD_LAZY_CALL_REF(Arity, Class)
* TRIQS_LAZY_ADD_LAZY_CALL_WITH_COPY(Arity, Class)

In both cases, calling the object will return a lazy expression. This operator() is 
enabled iif at least one argument is a lazy expression.

The two macro differ in the copy behaviour :

* TRIQS_LAZY_ADD_LAZY_CALL_REF : 
   
  The object is captured by reference.

* TRIQS_LAZY_ADD_LAZY_CALL_WITH_COPY

  The object is captured by value, i.e. a copy is made (or a view is made is some cases, Cf...).


Example ::

  struct MyObject { 

  // the non-lazy call operator
  double operator()(double x) const { return 10*x;}

  // Add call for lazy expression
  TRIQS_LAZY_ADD_LAZY_CALL_REF(1,MyObject);

  };

This object is then used as ::

  placeholder <2> x_;
  MyObject ob;
  auto ex =  2* ob(x_) + 3*x_; // ex is a lazy expression
  std::cout << ex<< std::endl; // prints the expression tree
  ex(4) ; // evaluates

* NB : when the non-lazy call operator is already a template

  The operator() created by the macro is enabled iif at least *one* argument is a lazy expression.
  If the other operator() are also template, they must be *disabled* in this case,
  using the `boost::disable_if` and the `triqs::clef::is_lazy` metafunctions. 
  Example::

    struct MyObject { 
     
    // This template is disabled is one argument is a lazy expression.
    template <typename Fnt>
    typename boost::disable_if< triqs::clef::one_is_lazy <Fnt>, double >::type 
    operator() (Fnt const & f ) const;
     
    // Add call for lazy expression
    TRIQS_LAZY_ADD_LAZY_CALL_REF(1,MyObject);
    };

  Or with two variables::

    struct MyObject { 

    // This template is disabled is one argument is a lazy expression.
    template <typename F1, typename F2>
    typename boost::disable_if< triqs::clef::one_is_lazy <F1,F2>, double >::type 
    operator() (F1 const & f1, F2 const & f2) const;
       
    // Add call for lazy expression
    TRIQS_LAZY_ADD_LAZY_CALL_REF(2,MyObject);
    };


