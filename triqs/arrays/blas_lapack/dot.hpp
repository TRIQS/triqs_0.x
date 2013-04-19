/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012 by O. Parcollet
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
#ifndef TRIQS_ARRAYS_BLAS_LAPACK_DOT_H
#define TRIQS_ARRAYS_BLAS_LAPACK_DOT_H
#include "./tools.hpp"
//#include "./qcache.hpp"

namespace triqs { namespace arrays { namespace blas {

 namespace f77 { // overload

  extern "C" { 
   double TRIQS_FORTRAN_MANGLING(ddot)(const int & N, const double * x, const int& incx, const double * y, const int & incy);
   double TRIQS_FORTRAN_MANGLING(zdotc)(const int & N, const std::complex<double> * x, const int& incx, const std::complex<double> * y, const int& incy);
   double TRIQS_FORTRAN_MANGLING(zdotu)(const int & N, const std::complex<double> * x, const int& incx, const std::complex<double> * y, const int& incy);
  }

  inline double dot (const int & M, const double* x, const int & incx, const double* Y, const int & incy)  { 
   return TRIQS_FORTRAN_MANGLING(ddot)(M, x, incx, Y, incy);
  }
  //inline std::complex<double> dot (const int & M, const std::complex<double>* x, const int & incx, const std::complex<double>* Y, const int & incy)  {
  // return TRIQS_FORTRAN_MANGLING(zdotu)(M, x, incx, Y, incy);
  // }

 }

 /**
  * Calls dot product of 2 vectors.
  * Takes care of making temporary copies if necessary
  */
 template<bool Star, typename VTX, typename VTY>
  typename std::enable_if< std::is_same<typename VTX::value_type,double>::value && have_same_value_type< VTX, VTY>::value, typename VTX::value_type >::type
  dot (VTX const & X, VTY const & Y) { 
   static_assert( is_amv_value_or_view_class<VTX>::value, "blas1 bindings only take vector and vector_view");
   static_assert( is_amv_value_or_view_class<VTY>::value, "blas1 bindings only take vector and vector_view");
   if (( X.size() != Y.size()) ) TRIQS_RUNTIME_ERROR << "Dimension mismatch in dot : X : "<<X().shape()<<" and Y : "<<Y().shape();
   //const_qcache<VTX> Cx(X); // mettre la condition a la main
   //const_qcache<VTY> Cy(Y); // mettre la condition a la main
   return f77::dot(X.size(), X.data_start(), X.stride(), Y.data_start(), Y.stride());
   //return f77::dot(X.size(), Cx().data_start(), Cx().stride(), Cy().data_start(), Cy().stride());
  }

 template< bool Star, typename T>
  typename std::enable_if<boost::is_complex<T>::value && Star,T>::type
  _conj(T && x) { return conj(std::forward<T>(x));}

 template< bool Star, typename T>
  typename std::enable_if<!( boost::is_complex<T>::value && Star),T>::type
  _conj(T && x) { return std::forward<T>(x);}

 /**
  * Calls dot product of 2 vectors.
  * Takes care of making temporary copies if necessary
  * general case. Also for complex since there is a bug on some machines (os X, weiss...) for zdotu, zdotc...
  * a transcription from netlib zdotu
  */
 template< bool Star, typename VTX, typename VTY>
  typename std::enable_if< (!std::is_same<typename VTX::value_type,double>::value && have_same_value_type< VTX, VTY>::value),
	   decltype(std::declval<VTX>()(0)* std::declval<VTY>()(0)) >::type
	   //decltype(std::declval<typename VTX::value_type>()* std::declval<typename VTY::value_type>()) >::type
  dot (VTX const & X, VTY const & Y) {
   if (( X.size() != Y.size()) ) TRIQS_RUNTIME_ERROR << "Dimension mismatch in dot : X : "<<X().shape()<<" and Y : "<<Y().shape();
   //const_qcache<VTX> Cx(X); // mettre la condition a la main
   //const_qcache<VTY> Cy(Y); // mettre la condition a la main
   size_t N= X.size(), incx = X.stride(), incy = Y.stride();
   decltype(X(0)*Y(0)) res = 0;
   // This only works for object with data (ISP), not only from the concept...
   auto * restrict X_ = X.data_start();
   auto * restrict Y_ = Y.data_start();
   if ((incx==1) && (incy==1)) {
    for (size_t i=0; i<N; ++i) res += _conj<Star>(X_[i]) * Y_[i];
   }
   else { // code for unequal increments or equal increments  not equal to 1
    for (size_t i=0, ix=0, iy=0; i<N; ++i, ix += incx, iy +=incy) {res += _conj<Star>(X_[ix]) * Y_[iy]; }
   }
   // general code for the concept. Is it really slower ?
   //for (size_t i=0; i<N; ++i) res += _conj<Star>(X(i)) * Y(i);
   return res;
  }
}

template <typename VTX,typename VTY> auto dot  (VTX const & X, VTY const & Y) DECL_AND_RETURN( blas::dot<false>(X,Y));
template <typename VTX,typename VTY> auto dotc (VTX const & X, VTY const & Y) DECL_AND_RETURN( blas::dot<true>(X,Y));

}}// namespace

#endif

