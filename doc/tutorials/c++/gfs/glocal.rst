Local Green's functions
########################

Declaring a Green's function
-------------------------------
.. compileblock:: 

 
    #include <triqs/gf/imfreq.hpp>
    using triqs::gf;
    int main(){
     double beta=10;
     int Nfreq =100;
     auto fermionic_imfreq_mesh = make_gf_mesh<imfreq>(beta,triqs::gf::Fermion,Nfreq); 
     auto G_weiss_up = make_gf<imfreq>(fermionic_imfreq_mesh, make_shape(1,1), tail_t(1,1));
    }
