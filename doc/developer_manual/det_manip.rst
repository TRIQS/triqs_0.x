det_manip
===============

Introduction
--------------

The purpose of this little class is to regroup  standard block manipulations on determinant, used in several 
QMC.

Given a function :math:`F(x,y)`, and two sets of values of x and y :math:`x_i,y_i \ 0\leq i \leq N`,
we can construct the square :math:`N\times N` matrix 

.. math:: 
   
   M_{i,j} = F(x_i,y_j)

The goal of this class is to quickly compute the determinant and the inverse of this matrix.

More precisely, the class contains : 

* Data : 

  * a vector containing  :math:`x_i,y_i \ 0\leq i \leq N`
  * :math:`M^{-1}` and :math:`det M`

* Methods to quickly update the det and inverse of the matrix when one :  

  * adds/removes a line and a column (i.e. adding or removing one x and one y)
  * changes a line/colum, etc... 


Parameter & construction
-----------------------------

The template parameter is the **FunctionType**, the type of F,
which is completely general but has to model the concept (NECESSARY ?)

  * return_type : type returned by the function 
  * argument_type : type of the argument of the function

Operations
-----------------------------

* The possible operations on the matrix M are : 

  +------------+------------------------------+
  | Operation  | Effect                       |
  +============+==============================+
  | insert     | adding a line and a column   |
  +------------+------------------------------+
  | remove     | removing a line and a column |
  +------------+------------------------------+
  | change_col | changing one *y*             |
  +------------+------------------------------+
  | change_raw | changing one *x*             |
  +------------+------------------------------+


* Each operation *OP* is called in two steps: 

  * value_type try_OP(arguments ...) 

    Returns the ratio 

    .. math:: \frac{det M'}{det M}

    where M' would be the matrix after the operation is completed.

    try_OP **does NOT** modify the matrix :math:`M`.

  * void complete_operation() 

    Complete the last operation OP (the last called try_OP), by updating the list of x and y 
    and the inverse of the matrix to :math:`(M')^{-1}`.

* This structure is designed to write  Monte Carlo algorithms : 
  
  * the try part of the move calls some try_OP
  * if and only if the move is accepted, is the complete_operation called.

Under the hood ...
-------------------------

* All matrix algebra is made with BLAS calls.

* The storage is done in a compact way: when a column or row is added, 
  no data are shifted, it is added at the end of the matrix.
  However, the permutation of row and columns are handled by this class
  so that this is transparent for the user.

Example 
-------------





