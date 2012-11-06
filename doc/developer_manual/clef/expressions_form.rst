
.. highlight:: c

Forming lazy expressions
===========================


The first step consists in forming lazy expressions, and then evaluate them.
Expressions are made of :

* Placeholders
* Binary operations on expressions `(+, -, *, /, >, <, >=, <=, ==)`
* Callable objects (see below) called on expressions
* Conditional if_else expressions
 
Placeholders
----------------

Loosely speaking, a placeholder is a "variable name" used to build an expression.
Technically, it is a trivial object, with a type but containing no data.
Placeholders are declared as ::

  placeholder<Number> Name;

This declares a placeholder called Name (an empty object for C++). 

Example ::

  placeholder <1> x_; 
  placeholder <2> y_; 

Note that the only thing of significance in a placeholder is its type, since it contains no data.
In other words, a placeholder is **empty**. It contains **no value** at runtime. 
   
  .. warning:: 
    
      As a consequence, defining ::
      
        placeholder <1> y_; 

      would imply that `x_` is the same as `y_` : `x_` == `y_` will be always true.

Usual operators
----------------------

Standard arithmetic operations are possible with expressions, e.g.::
 
  auto e1 = x_ + 2* y_;

Callable objects
--------------------

* Objects can overload the operator () for lazy expressions in order to build more complex
  expressions.

* For example, the header `math.hpp` contains the declaration to make 
  the basic function of std `math.h` accept lazy_expressions.
  
 .. compileblock::
 
   #include <triqs/clef/math.hpp>
   int main () { 
    triqs::clef::placeholder <1> x_; 

    auto e1 = cos(2*x_+1);
    auto e2 = abs(2*x_-1);
    auto e3 = floor(2*x_-1);
    auto e4 = pow(2*x_+1,2);
   }
* To make your object callable, or to overload a function to accept lazy argument,  see :ref:`callable_object`.

*lazy* function
-------------------

The *lazy* function can be called on any object to make it lazy, e.g. 

.. compileblock::
 
   #include <triqs/clef/core.hpp>
   #include <vector>
   namespace tql = triqs::clef;
   int main () { 
    std::vector<int> V;
    tql::placeholder<1> i_;
    auto e1 = tql::lazy(V)[i_];
   }

NB : The type of expressions
-------------------------------

The type of expressions (e1,e2, ... above) is quite complex. The structure of the expression is actually encoded in the type
itself, using recursive template (a technique called *expression template*).
Therefore it worth noting that : 

* We use auto to declare an expression (*NB : non C++11 users should use BOOST_AUTO*) ::
  
   auto e1 = x_ + 2* y_;

* Declaring an expression does not do any computation, hence the name of the library (lazy).
  It just stores the expression tree (its structure in the type, and the leaves of the tree).


Copy policy in building expressions
---------------------------------------------------

A central question when forming expressions is whether the object at the leaves of the expressions tree
(scalar, placeholders, various callable objects, etc...) should be captured by value or by reference.

In the lazy library, the choice has been made to capture them **by value**, i.e. : 

  *By default, all objects appearing in a lazy expression are* **copied**, *rather than captured by reference*.

This is necessary to store expressions (with auto like above) for future reuse, transform them into new expressions
(e.g. make partial evaluation). Expressions are ordinary objects. 
If the leaves of the expression tree were captured by reference, a guarantee would have to be made that 
they will live at least as long as the expression itself, or one gets dangling references.

The drawback of this approach is that it can generate unless copies of large objects.
There are several solutions to this issue : 

* If the object is very small (like a double), hence making a copy in not a problem.

* If you *know* that the object `A` will survive the expression, so using a reference is not a problem.
  You can use the `lazy(A)` expression that will wrap the reference to the object.

* If the object has a compagnon view object (like array, array_view). In this case, 
  one wishes to put a view of the object rather than a copy in the expression.
  There are two sub-cases : 

  * TRIQS objects like array, matrix, vector will do it automatically.

  * For another object, if the object defines the tag `has_view_type_tag` as ::

     typedef void has_view_type_tag;
     typedef MY_VIEW_TYPE view_type;
  
    the instance of the object will be replaced by an instance of its view_type in building the expression.

For an illustration, Cf....


