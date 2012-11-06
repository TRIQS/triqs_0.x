
.. highlight:: c

Evaluating lazy expressions
===============================

Forming expressions is nice, but completely useless unless one can *evaluate* them 
by affecting some values to the placeholder and actually *do* the computation.

Expressions are evaluated with the eval function ::

 eval ( expression , placeholder_1 = value_1, placeholder_2 = value_2, ...)

The evaluation can be :

* *complete* : when a value has been affected to every placeholder, so the result is not
  an expression any more, but the result of the computation....

* *partial* : when some placeholder are still free (no value has been affected to them).
  The result is another expression where the parts that could be computed have been computed.

Complete evaluation
--------------------

Example ::

  eval (x_  + 2*y_ , x_=1, y_ = 2)

  auto e = x_  + 2*y_;
  auto r = eval( e , x_=1, y_ = 2);

Note that : 

* The order of placeholder does not matter in calling eval.
* It is an error to put the same placeholder twice.

.. note::
   gcc 4.4 seems to have an ADL bug. So you need to qualify eval as triqs::clef::eval in that case...

Partial evaluation
--------------------

The evaluation can also be partial, in which case the compiler replaces the placeholder whose value is provided
and rebuild a new expression tree.

Example ::

 auto e1 = x_  + 2*y_;
                                          // return type  ---> value 
 auto r  = eval( e1 , x_ =1);             // expression   ---> 1 + (2 * y_)
 auto e2 = eval( e1 , x_ =x_ + y_);       // expression   ---> ((x_ + y_) + (2 * y_))
 auto r  = eval( e2,  x_ = 1, y_ = 2);    // int          ---> 7

During the partial evaluation, all subtrees that can evaluated are evaluated.
For example ::

 eval(x_ + f(y_) , y_=1);                 // --> x + f(1)

In this case, the function f is called by eval.


