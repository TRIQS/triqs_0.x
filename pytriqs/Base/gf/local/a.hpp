#include <triqs/gf/matsubara_freq.hpp> 
#include <triqs/gf/matsubara_time.hpp> 
#include <triqs/gf/block.hpp> 

namespace tql= triqs::clef;
namespace tqa= triqs::arrays;
using tqa::range;
using triqs::arrays::make_shape;
using triqs::gf::gf;
using triqs::gf::gf_view;
using triqs::gf::block;
using triqs::gf::Fermion;
using triqs::gf::matsubara_freq;
using triqs::gf::matsubara_time;

void test_block ( gf_view<block<matsubara_freq>> const & G) { 

 auto  g0v = G(0)();
 g0v(0) = 3.2;
 
}

void test_gf ( gf_view<matsubara_freq> const & G) { 

 auto  g0v = G();
 g0v(0) = -3.2;
 
}

void test2 ( tqa::array<double,2> const & M) { 

 std::cout  << "M = "<< M<< std::endl;
 
}

void test3 ( std::vector<tqa::matrix_view<double> > & M) { 

 std::cout  << "M = "<< M[0]<< std::endl;
 
 M[0] (0,0) = -100;
}

