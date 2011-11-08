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

#include "./src/matrix.hpp"
#include "./src/linalg/matmul.hpp"
#include "./src/linalg/mat_vec_mul.hpp"
#include "./src/linalg/inverse.hpp"
#include "./src/linalg/determinant.hpp"

using namespace triqs::arrays;
const int N1= 200, N2 = 300;

struct plain_for_no_ptr { 
 void operator()() { 
  triqs::arrays::matrix<double,Option::Fortran > A (N1,N2);
  for (int u =0; u<1000; ++u)
   for (int i =0; i<N1; ++i)
    for (int j=0; j<N2; ++j) 
     A(i,j) = 10*i+ j;
 }
};

template<class int_type>
struct pointer_generic_form { 
 void operator()() { 
  triqs::arrays::matrix<double,Option::Fortran > A (N1,N2);
   const int_type s0 = A.indexmap().strides()[0];
   const int_type s1 = A.indexmap().strides()[1];

   const int_type l0 = A.indexmap().lengths()[0];
   const int_type l1 = A.indexmap().lengths()[1];
  for (int u =0; u<1000; ++u) { 
   double * restrict p= &(A(0,0));
   //typedef int int_type; 

   //assert(s0==1); assert(s1==N1);
   //assert(l0==N1); assert(l1==N2);
   for (int_type i0 =0; i0<l0; ++i0)
    for (int_type i1 =0; i1<l1; ++i1)
     p[i0*s0 + i1*s1] = 10*i0+ i1;

  }
 }
};

struct iterators { 
 void operator()() { 
  typedef triqs::arrays::matrix<double,Option::Fortran > MM; 
  MM A (N1,N2);
  for (int u =0; u<1000; ++u)
  {
   for (MM::iterator it = A.begin(); it; ++it) { 
    *it =10*it.indices()[0] + it.indices()[1] ;
   }
  }
 }
};

#include "./src/indexmaps/cuboid/foreach.hpp"

typedef double value_type;
typedef triqs::arrays::matrix<double,Option::Fortran >::indexmap_type::domain_type::index_value_type index_value_type;
struct F { 
void operator()(value_type & p, index_value_type const & key) const { p=key[0]*10 + key[1] ;}
};

struct foreach  {
 void operator()() {
  triqs::arrays::matrix<double,Option::Fortran > A (N1,N2);
  triqs::arrays::indexmaps::foreach(F(),A);
 // for (int u =0; u<5000; ++u)    make_view(A) = 1867;
 }
};


struct const_with_iterators {
 void operator()() {
  typedef triqs::arrays::matrix<double,Option::Fortran > MM;
  MM A (N1,N2);
  for (int u =0; u<5000; ++u)
  {
   for (MM::iterator it = A.begin(); it; ++it) {
    *it = 1876;
   }
  }
 }
};

struct plain_for_no_ptr_const { 
 void operator()() { 
  triqs::arrays::matrix<double,Option::Fortran > A (N1,N2);
  for (int u =0; u<5000; ++u)
   for (int i =0; i<N1; ++i)
    for (int j=0; j<N2; ++j) 
     A(i,j) = 1876;
 }
};

struct assign_to_const { 
 void operator()() { 
  triqs::arrays::matrix<double,Option::Fortran > A (N1,N2);
  for (int u =0; u<5000; ++u)
     make_view(A) = 1867;
 }
};



#include "./speed_tester.hpp"
int main() {
 const int l = 100;
 speed_tester<pointer_generic_form <std::ptrdiff_t> > (l);
 speed_tester<pointer_generic_form < size_t> > (l);
 speed_tester<plain_for_no_ptr> (l);
 speed_tester<foreach> (l);
 speed_tester<plain_for_no_ptr_const> (l);
 speed_tester<assign_to_const> (l);
 speed_tester<iterators> (l);
 speed_tester<const_with_iterators> (l);
 return 0;
}

