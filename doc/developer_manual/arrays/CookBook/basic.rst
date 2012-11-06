Basic 
============

.. highlight:: c


Declaring and printing an array
-------------------------------
.. compileblock:: 

 
    #include <triqs/arrays.hpp>
    namespace tqa=triqs::arrays;
    int main(){
       tqa::array<double,1> A(20);
       std::cout << "A = "<<A << std::endl; // arrays are not init by default: this is random 
       A() = 2;     //  assign 2 to a complete view of A.
       std::cout <<"A = "<< A << std::endl;
       A(4) = 5;
       std::cout <<"A = "<< A << std::endl;

    }


Simple operations
-------------------

.. compileblock:: 

    #include <triqs/arrays.hpp>
    namespace tqa=triqs::arrays;
    int main(){
      tqa::array<double,1> A(10),B(10);
      A()=2;B()=3;
      tqa::array<double,1> C = A+B;
      std::cout << "C = "<<C << std::endl;
    }


Simple functions
-------------------
.. compileblock:: 

    #include <triqs/arrays.hpp>
    namespace tqa=triqs::arrays;
    int main(){
      tqa::array<double,1> A(10),B(10);
      A()=2;B()=3;
      ///.....
    }


HDF5 Archiving
-------------------
Archiving an array into an HDF5 file is easy:

.. compileblock::

    #include <triqs/arrays.hpp>
    namespace tqa=triqs::arrays;
    int main(){
    
      tqa::array<double,2> A(2,2); A() = 3;          // declare and init

      H5::H5File file("store_A.h5",H5F_ACC_TRUNC);   // open the file
      h5_write(file,"A",A);                         // write the array as 'A' into the file

      tqa::array<double,2> B;                        // read the file into B
      h5_read (file, "A",B);               
      std::cout << "B = "<<B<<std::endl;
    }


Views: ranges and slices
-------------------------
One can easily take a slice of an array to view and modify only part of the underlying data.

.. compileblock::

    #include <triqs/arrays/array.hpp>
    namespace tqa=triqs::arrays;
    int main(){
      tqa::array<double,2> A(3,3); A() = 2.5;   
      std::cout << A <<std::endl;
      
      tqa::array_view<double,1> B = A(1,tqa::range()); //select the first line of the matrix
      std::cout <<"B = "<< B << std::endl;
      B(0) = 1;

      std::cout <<"A = "<< A << std::endl;            
    }


Matrices and vectors
-------------------------
Arrays must be distinguished from vectors and matrices, which have an algebra of their own.

.. compileblock::
    
    #include <triqs/arrays.hpp>

    namespace tqa=triqs::arrays;
    int main(){
     tqa::array<double,2> A(2,2), B(2,2),C; 
     
     A() = 3; B() = 1; C = A*B;
     std::cout << "A*B = "<< C << std::endl;

     tqa::matrix<double> D(2,2),E(2,2),F; 
     E() = 3; E() = 1; F = D*E;
     std::cout << "C*D = "<< F << std::endl;

     tqa::vector<double> u(2),v(2),w;
     u()=1;v()=2; w = u+v;
     
     std::cout <<"u+v = "<< w << std::endl;
    }



Defining through a lazy expression
-----------------------------------

.. compileblock::

    #include <triqs/arrays.hpp>
    namespace tqa=triqs::arrays; namespace tql=triqs::clef;
 
    int main(){
       tql::placeholder<0> i_;   tql::placeholder<1> j_;
       tqa::array<double,2> A(2,2);  
       A(i_,j_)= i_ + j_ ;
       std::cout << "A = "<<A << std::endl;
    }



Linear algebra
---------------

.. compileblock::

    #include <triqs/arrays.hpp>
    #include <triqs/arrays/linalg/inverse.hpp>
    #include <triqs/arrays/linalg/determinant.hpp>
    
    namespace tql=triqs::clef; namespace tqa=triqs::arrays;
    int main(){
      tql::placeholder<0> i_;
      tql::placeholder<1> j_;
      tqa::matrix<double> A(2,2); 

      A(i_,j_) = i_+j_; 
      tqa::matrix<double> B = inverse(A); 
      double C = determinant(A); 
 
      std::cout << "A^(-1) = "<< B << std::endl;
      std::cout << "det(A) = " <<C <<std::endl;

    }


Map and fold
-------------

.. compileblock::
  
    #include <triqs/arrays.hpp>
    #include <triqs/arrays/functional/map.hpp>
    namespace tqa=triqs::arrays;
    
    double f(int i) { return i*10;}

    int main() {
      auto F = tqa::map(boost::function<double(int)>(f));
      tqa::array<int,2> A(2,2); A() =2;
 
      tqa::array<double,2> B,C;

      A() =2;
      B = F(A);
      C = F( 2*A );  // works also with expressions of course


      std::cout << "A = "<<A<<std::endl;
      std::cout << "F(A) = "<<B<<std::endl;
      std::cout << "F(2*A) = "<<C<<std::endl;

    }


Bound checking
---------------
By default, there is no bound checking:

.. compileblock::

    #include <triqs/arrays.hpp>
    namespace tqa=triqs::arrays;
    int main(){
        tqa::array<double,2> A(2,2); A() = 3;   
        std::cout << A(0,3) << std::endl;            
    }

But one can add bound-checking by adding a preprocessor command:

.. compileblock::

    #define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK
    #include <triqs/arrays.hpp>
    namespace tqa=triqs::arrays;
    int main(){
        tqa::array<double,2> A(2,2); A() = 3;   
        std::cout << A(0,3) << std::endl;            
    }



  
