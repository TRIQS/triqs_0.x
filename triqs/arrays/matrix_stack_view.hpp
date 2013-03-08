/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013 by O. Parcollet
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
#ifndef TRIQS_ARRAYS_MATRIX_STACK_VIEW_H
#define TRIQS_ARRAYS_MATRIX_STACK_VIEW_H
#include "./array.hpp"
#include "./matrix.hpp"
#include "./matrix_view_proxy.hpp"
#include <triqs/arrays/linalg/inverse.hpp>

namespace triqs { namespace arrays {

 template<typename T> class matrix_stack_view { 
  public:
   typedef array_view<T,3> array_view_t;

   matrix_stack_view( typename array_view_t::view_type const & a_):a(a_) {}
   matrix_stack_view( array_view_t && a_):a(std::move(a_)) {}
   
   matrix_view_proxy      <array_view_t,2> operator()(size_t i)       { return matrix_view_proxy       <array_view_t,2>(a,i);}
   const_matrix_view_proxy<array_view_t,2> operator()(size_t i) const { return const_matrix_view_proxy <array_view_t,2>(a,i);}
 
   matrix_view<T> view(size_t i) const { return a(range(),range(),i);}
   
   size_t dim0() const { return a.len(0);}
   size_t dim1() const { return a.len(1);}
   size_t size() const { return a.len(2);}

   matrix_stack_view & operator +=(matrix_stack_view const & arg) { a += arg.a; return *this; }
   matrix_stack_view & operator -=(matrix_stack_view const & arg) { a -= arg.a; return *this; }

   template<typename RHS>   
    typename std::enable_if<RHS::rank ==2, matrix_stack_view &>::type 
    operator +=(RHS const & arg) { for (size_t i=0; i<size(); ++i) view(i) +=arg; return *this; }

   template<typename RHS>   
    typename std::enable_if<RHS::rank ==2, matrix_stack_view &>::type 
    operator -=(RHS const & arg) { for (size_t i=0; i<size(); ++i) view(i) -=arg; return *this; }
  
   template<typename RHS> matrix_stack_view & operator *=(RHS const & arg) { a*= arg; return *this; }
   template<typename RHS> matrix_stack_view & operator /=(RHS const & arg) { a/= arg; return *this; }

   void invert() {for (size_t i=0; i<size(); ++i) { auto v = view(i); v = inverse(v);} }

   friend matrix_stack_view matmul_R_L ( matrix_view<T> const & L, matrix_stack_view const & M, matrix_view<T> const & R) { 
    matrix_stack_view res (typename array_view_t::non_view_type (L.dim0(), R.dim1(),M.size()));
    for (size_t i=0; i<M.size(); ++i)  { res.view(i) = L * M.view(i) * R; }
    return res;
   }

  private:
   array_view_t a;
 };

}}
#endif

