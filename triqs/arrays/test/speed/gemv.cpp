/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011 by O. Parcollet
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#include <boost/numeric/bindings/std/vector.hpp>
#include <boost/numeric/bindings/blas/level1/axpy.hpp>
#include <boost/numeric/bindings/blas/level2/ger.hpp>
#include <boost/numeric/bindings/blas/level2/gemv.hpp>
#include <boost/numeric/bindings/blas/level3/gemm.hpp>
#include <boost/numeric/bindings/lapack/computational/getrf.hpp>
#include <boost/numeric/bindings/lapack/computational/getri.hpp>

#include "./src/vector.hpp"
#include "./src/matrix.hpp"
#include "./src/proto/matrix_algebra.hpp"
#include "./src/linalg/matmul.hpp"
#include "./src/linalg/mat_vec_mul.hpp"

using namespace std;
using namespace triqs;
using namespace triqs::arrays;
//using namespace triqs::arrays::linalg;
namespace blas = boost::numeric::bindings::blas;
namespace lapack = boost::numeric::bindings::lapack;
namespace bindings= boost::numeric::bindings;

struct gemv_via_binder { 

 triqs::arrays::matrix<double> A;
 typedef triqs::arrays::vector<double> vector_type;
 vector_type  MC, MB;
 static const unsigned int N =200;

 gemv_via_binder(): A(N,N,FORTRAN_LAYOUT), MC(N), MB(N) {
 
  for (int i =0; i<N; ++i)
    for (int j=0; j<N; ++j) 
     A(i,j) = 0.1*(10*i+ j);

   make_view(MC) = 1;

 }

 void operator()() { 
  blas::gemv(1,A, MC, 0, MB);
 }
};

struct expressif : gemv_via_binder  { 
 void operator()() { 
  MB = A * MC;
 }
};

namespace myblas { 
//#include "./myblas.hpp"
}

#define MYFORTRAN( id ) id##_

//void MYFORTRAN(dgemv)(const char* trans, const int & m, const int & n, const double & alpha, const double A[], int & lda,
//			const double x[], const int & incx, const double & beta, double y[], const int & incy);
 
struct plain_gemv_call : gemv_via_binder  { 
 void operator()() { 

 const char CT ='N';
 const double Done=1,Dzero=0;
 const int one =1;
 double * A_ = &A(0,0);
 double * MB_ = &MB(0);
 double * MC_ = &MC(0);

 int n = N;
 //myblas::MYFORTRAN(dgemv)(&CT,n,n,Done, A_,n,MC_ ,one,Dzero,MB_,one);
 MYFORTRAN(dgemv)(&CT,&n,&n,&Done, A_,&n,MC_ ,&one,&Dzero,MB_,&one);
 }
};


#include "./speed_tester.hpp"
int main() {
 const int l= 1000*1000;
 speed_tester<plain_gemv_call> (l);
 speed_tester<gemv_via_binder> (l);
 speed_tester<expressif> (l);
 return 0;
}

