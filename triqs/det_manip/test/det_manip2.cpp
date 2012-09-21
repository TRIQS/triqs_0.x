#include <iostream>
#include <triqs/det_manip/det_manip.hpp>
#include <triqs/mc_tools/random_generator.hpp>
#include <triqs/arrays/linalg/det_and_inverse.hpp>
#include <triqs/arrays/impl/asserts.hpp>

struct fun {
 typedef double result_type;
 typedef double argument_type;
 double operator()(double x, double y) const { return 0.5* ( y>=x ? 1 : -1) ; }
};

template<class T1, class T2 >
void assert_close( T1 const & A, T2 const & B, double precision) {
 if ( abs(A-B) > precision) TRIQS_RUNTIME_ERROR<<"assert_close error : "<<A<<"\n"<<B;
}
const double PRECISION = 1.e-12;

struct test { 

 fun f;
 triqs::det_manip::det_manip<fun> D;
 double det_old,detratio;

 test() : f(), D(f,100) {}

 void check() { 
#ifndef PRINT_ALL
  std::cerr  << "det = " << D.determinant() <<  " == " << double(determinant(D.matrix()))<< std::endl;
#else
  std::cerr  << "det = " << D.determinant() <<  " == " << double(determinant(D.matrix()))<< std::endl <<
   D.inverse_matrix() << D.matrix() << triqs::arrays::matrix<double>(inverse(D.matrix()))<<  std::endl;
#endif
  assert_close(D.determinant() , 0.5, PRECISION); 
  assert_close(D.determinant() , 1/determinant(D.inverse_matrix()), PRECISION); 
  triqs::arrays::assert_all_close( inverse(D.matrix()) , D.inverse_matrix(), PRECISION);
  assert_close( det_old * detratio , D.determinant(), PRECISION);
 }

 void run() { 
  triqs::mc_tools::random_generator RNG("mt19937", 23432);
  for (size_t i =0; i< 100; ++i) { 
   std::cerr <<" i = "<< i << " size = "<< D.size() << std::endl; 
   // choose a move
   int mn = RNG(4);
   size_t s = D.size();
   size_t w;
   det_old = D.determinant();
   double x,y; 

   switch(RNG(( i>10 ? 2 : 1))) {
    case 0 :
     x = RNG(10.0), y = RNG(10.0);
     w = RNG(s);
     detratio = D.try_insert(w,w, x,x); 
     break;
    case 1 : 
     w = RNG(s);
     if (s>0) detratio = D.try_remove(w,w);
     break;
    default : 
     TRIQS_RUNTIME_ERROR <<" TEST INTERNAL ERROR" ;
   };

   D.complete_operation();
   if (D.size() >0) check();
  }
 }

};

int main(int argc, char **argv) {
 test().run();
}


