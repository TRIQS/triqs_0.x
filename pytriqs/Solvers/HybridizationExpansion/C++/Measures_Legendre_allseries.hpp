
/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011 by M. Ferrero, O. Parcollet
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

#ifndef TRIQS_CTHYB1_MEASURES_LEGENDRE_ALL_H
#define TRIQS_CTHYB1_MEASURES_LEGENDRE_ALL_H
#include <triqs/arrays/array.hpp>
#include <triqs/arrays/h5/array_stack.hpp>
#include <triqs/gf_local/GF_Bloc_ImLegendre.hpp>
#include "Measures_Z.hpp"
#include "Configuration.hpp"

namespace tqa=triqs::arrays;

/* Measure in Legendre Polynoms.  */
class Measure_G_Legendre_all : public Measure_acc_sign<COMPLEX> { 

 const Configuration & conf;
 const int a_level;
 GF_Bloc_ImLegendre & Gl;
 typedef Measure_acc_sign<COMPLEX> BaseType;  

 tqa::array<double,3> data,data_part;
 tqa::h5::H5File outfile;
 tqa::h5::array_stack< tqa::array<double,3 > >  data_stack;

 template<class T1, class T2, class T3, class T4>
  static std::string filename(T1 x1, T2 x2, T3 x3, T4 x4) { std::stringstream f; f<<x1<<x2<<x3<<x4; return f.str(); }

 public:
 const std::string name;  

 // PUT the file out to keep the node...
 // multinode....
 Measure_G_Legendre_all (Configuration const & conf_, int a, GF_Bloc_ImLegendre & Gl_):
  BaseType(), conf(conf_), a_level(a), Gl(Gl_), name(to_string("G(legendre)",a)), 
 //data (Gl.data.shape()), // need triqs::arrays for this
 //data_part (Gl.data.shape())
  data(Gl.N1, Gl.N2, Gl.mesh.index_max+1),  
  data_part(Gl.N1, Gl.N2, Gl.mesh.index_max+1) ,  
  outfile(filename("Gl_stack",a,".h5",""), H5F_ACC_TRUNC ) ,
  data_stack(outfile, "Gl",data.shape(), 1000)
 {
  assert( Gl.mesh.index_min ==0);
  data()=0;
 }

 void accumulate(COMPLEX signe);

 void collect_results( boost::mpi::communicator const & c){
  //std::cerr<< "finalizing "<< name <<std::endl;
  //std::cerr<< " h5stack size "<< data_stack.size() <<std::endl;
  data_stack.flush();
  outfile.close(); 
  BaseType::collect_results(c);
  double Z_qmc ( real(this->acc_sign));

  return ;
  // not to collide with the old measure
  for (int n1=1; n1<=Gl.N1;n1++)
   for (int n2=1; n2<=Gl.N2;n2++)
    for (int l = Gl.mesh.index_min; l <= Gl.mesh.index_max; l++)
     Gl.data(n1,n2,l) = -(sqrt(2.0*l+1.0)/(Z_qmc*conf.Beta)) * data(n1-1,n2-1,l);
  Gl.MPI_reduce_sum_onsite();
  Gl.MPI_bcast();
  Gl.determine_tail();

 }

};

#endif
