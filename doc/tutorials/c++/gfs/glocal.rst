Local Green's functions
########################

Declaring a Green's function
-------------------------------
.. compileblock:: 

 
    #include <triqs/gf/imfreq.hpp>
    using namespace triqs::gf;
    int main(){
     double beta=10;
     int Nfreq =100;
     auto fermionic_imfreq_mesh = make_gf_mesh<imfreq>(beta,triqs::gf::Fermion,Nfreq); 
     auto G_weiss_up = make_gf<imfreq>(fermionic_imfreq_mesh, triqs::arrays::make_shape(1,1), local::tail(1,1));
    }
