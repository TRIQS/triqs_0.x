#include <triqs/gf/local/gf.hpp>

//using namespace triqs::gf::local;
using namespace triqs::gf;
namespace tql= triqs::lazy;

int main() {

 std::vector<std::string> ind; ind.push_back("1");ind.push_back("2");
 local::gf<meshes::matsubara_freq,false> G(2,2, meshes::matsubara_freq(Fermion), ind,ind);
 //local::gf<meshes::matsubara_freq,true>  Ge(2,2, meshes::matsubara_freq(Fermion), ind,ind);

 //local::gf<meshes::matsubara_freq,false> Gd; 
 //local::gf<meshes::matsubara_freq,true> Gvd; 
 
 local::gf<meshes::matsubara_freq,true> Gv =G;
 std::cout  << G( 0) <<std::endl ;

 triqs::lazy::placeholder<0> om_;

 std::cout << G(om_) << std::endl ;
 std::cout << tql::eval(G(om_), om_=0) << std::endl ;

 std::cout << Gv(om_) << std::endl ;
 std::cout << tql::eval(Gv(om_), om_=0) << std::endl ;

 G(om_) = 1/(om_ + 2.3);

  std::cout << Gv(om_) << std::endl ;
 std::cout << tql::eval(Gv(om_), om_=0) << std::endl ;

 
}
