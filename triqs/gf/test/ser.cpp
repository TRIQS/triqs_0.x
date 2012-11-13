//#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK

#include <triqs/utility/serialization.hpp>
#include <triqs/gf/imfreq.hpp> 
#include <triqs/gf/imtime.hpp> 

namespace tql= triqs::clef;
namespace tqa= triqs::arrays;
using tqa::range;
using triqs::arrays::make_shape;
using triqs::gf::Fermion;
using triqs::gf::imfreq;
using triqs::gf::imtime;

#define TEST(X) std::cout << BOOST_PP_STRINGIZE((X)) << " ---> "<< (X) <<std::endl<<std::endl;

int main() {

 triqs::gf::freq_infty inf;

 triqs::gf::gf<triqs::gf::imfreq> G1; // empty
 TEST( G1( 0) ) ;

 double beta =1;
 auto G =  imfreq::make_gf (beta, Fermion, make_shape(2,2));

 double x = 127;
 std::string s = triqs::serialize(x);
 
 std::cout  << " s = "<< s<< std::endl;

 std::cout  << triqs::deserialize<double>(s) << std::endl;
 std::cout  << triqs::deserialize<int>(s) << std::endl;

std::vector<std::string> v; v.push_back("abc"); v.push_back("3");
std::cout << triqs::serialize(v)<< std::endl;

}
