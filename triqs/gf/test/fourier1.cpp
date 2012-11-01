//#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK

#include <triqs/gf/matsubara_freq.hpp> 
#include <triqs/gf/matsubara_time.hpp> 
#include <triqs/gf/local/fourier_matsubara.hpp> 

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

 double beta =1;
 auto G =  matsubara_freq::make_gf (beta, Fermion, make_shape(2,2));
 auto Gc = matsubara_freq::make_gf (beta, Fermion, make_shape(2,2));
 auto G3 = matsubara_freq::make_gf (beta, Fermion, make_shape(2,2));
 auto Gt = matsubara_time::make_gf (beta, Fermion, make_shape(2,2));

 auto gt = inverse_fourier(G);
 auto gw = fourier(gt);

 gw() = lazy_fourier(gt);
}


