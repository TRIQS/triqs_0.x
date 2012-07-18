#include "./python_stuff.hpp"
#include <triqs/clef/core.hpp>
#include <triqs/arrays/indexmaps/cuboid/group_indices.hpp>
#include <triqs/arrays/linalg/inverse.hpp>
#include <triqs/arrays/proto/array_algebra.hpp>
#include <triqs/clef/adapters/array.hpp>
#include <triqs/arrays/impl/asserts.hpp>

/*
   struct to_vec_ {
   std::ostream & out; to_vec_ (std::ostream & out_):out(out_){}
   template<typename T> void operator()(T) { out<< T(); }
   };

   template <typename V> void print() { boost::mpl::for_each<V>(to_vec_(std::cout)); }
   */
namespace tqa=triqs::arrays;
namespace tql=triqs::clef;
namespace mpl=boost::mpl;

using tqa::m_index;

template<typename STO> void test() { 

 tql::placeholder<0> i_;
 tql::placeholder<1> j_;
 tql::placeholder<2> k_;
 tql::placeholder<3> l_;

 { // a simple test
  tqa::array<int,4,STO> A(2,2,2,2); 
  A(i_,j_,k_,l_)= i_ + 10*j_ + 100*k_ + 1000*l_;
  TEST(A);
  TEST( group_indices(A,  m_index<0,1>(), m_index<2,3>())); 
 }

 { // more complex : inversing a tensor product of matrices...
  tqa::matrix<double,STO> B(2,2), C(3,3), Binv, Cinv;
  C(i_,j_) = 1.7 /( 3.4*i_ - 2.3*j_ + 1) ;
  B(i_,j_) = 2*i_ + j_;
  TEST(B); TEST(C);
  Binv = inverse(B);
  Cinv = inverse(C);
  TEST(Binv); TEST(Cinv);

  tqa::array<double,4,STO> A(2,3,2,3); 
  A(i_,j_,k_,l_) = B(i_,k_) * C(j_,l_);	

  TEST(A);

  tqa::matrix_view<double,STO> M = group_indices (A, m_index<0,1>() , m_index<2,3>() );
  M = inverse(M);

  // checking 
  tqa::array<double,4,STO> R(A.shape());
  R(i_,j_,k_,l_) = Binv(i_,k_) * Cinv(j_,l_);	

  TEST(R);
  TEST(A);
  assert_all_close(R,A,1.e-15);
  //TEST(make_array(R-A));
  //TEST( max(abs(R-A)));
 }

}

int main () {
 test<tqa::Option::Fortran>();
 test<tqa::Option::C>();

 // an alternative syntax ? Why to use the npl here ??
 //auto V = tqa::group_indices( A(i_,j_,k_,l_), m_index(i_,j_), m_index(k_,l_));

}
