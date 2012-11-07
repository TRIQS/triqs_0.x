//#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK

#include "../ser.hpp"
#include <triqs/gf/matsubara_freq.hpp> 
#include <triqs/gf/matsubara_time.hpp> 

namespace tql= triqs::clef;
namespace tqa= triqs::arrays;
using tqa::range;
using triqs::arrays::make_shape;
using triqs::gf::Fermion;
using triqs::gf::matsubara_freq;
using triqs::gf::matsubara_time;

#define TEST(X) std::cout << BOOST_PP_STRINGIZE((X)) << " ---> "<< (X) <<std::endl<<std::endl;

int main() {

 triqs::gf::freq_infty inf;

 triqs::gf::gf<triqs::gf::matsubara_freq> G1; // empty
 TEST( G1( 0) ) ;

 double beta =1;
 auto G =  matsubara_freq::make_gf (beta, Fermion, make_shape(2,2));

 double x = 127;
 std::string s = triqs::serialize(x);
 
 std::cout  << " s = "<< s<< std::endl;

 std::cout  << triqs::deserialize<double>(s) << std::endl;
 std::cout  << triqs::deserialize<int>(s) << std::endl;
}
