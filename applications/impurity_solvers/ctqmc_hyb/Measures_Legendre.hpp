
/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011 by L. Boehnke, M. Ferrero, O. Parcollet
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

#ifndef TRIQS_CTHYB1_MEASURES_LEGENDRE_H
#define TRIQS_CTHYB1_MEASURES_LEGENDRE_H
#include <triqs/gf_local/GF_Bloc_ImLegendre.hpp>
#include "Measures_Z.hpp"

/* Measure in Legendre Polynoms.  */
class Measure_G_Legendre : public Measure_acc_sign<COMPLEX> { 

 const Configuration & conf;
 const int a_level;
 GF_Bloc_ImLegendre & Gl;
 typedef Measure_acc_sign<COMPLEX> BaseType;  

 public:
 const std::string name;  

 Measure_G_Legendre (Configuration const & conf_, int a, GF_Bloc_ImLegendre & Gl_):
  BaseType(), conf(conf_), a_level(a), Gl(Gl_), name(to_string("G(legendre)",a)){Gl.data=0;}

 void accumulate(COMPLEX signe);

 void collect_results( boost::mpi::communicator const & c){
   BaseType::collect_results(c);
   double Z_qmc ( real(this->acc_sign));

   for (int n1=1; n1<=Gl.N1;n1++)
     for (int n2=1; n2<=Gl.N2;n2++)
       for (int l = Gl.mesh.index_min; l <= Gl.mesh.index_max; l++)
         Gl.data(n1,n2,l) = -(sqrt(2.0*l+1.0)/(Z_qmc*conf.Beta)) * Gl.data(n1,n2,l);
   Gl.MPI_reduce_sum_onsite();
   Gl.MPI_bcast();
   Gl.determine_tail();

 }

};

#endif
