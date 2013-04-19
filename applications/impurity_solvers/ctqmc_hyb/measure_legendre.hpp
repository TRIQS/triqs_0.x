
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

#include <triqs/gf/legendre.hpp>
#include "measure_z.hpp"

/* Measure in Legendre Polynoms.  */
class Measure_G_Legendre : public Measure_acc_sign<std::complex<double>> { 

 const Configuration & conf;
 const int a_level;
 triqs::gf::gf_view<triqs::gf::legendre> & Gl;
 typedef Measure_acc_sign<std::complex<double>> BaseType;  

 public:
 const std::string name;  

 Measure_G_Legendre (Configuration const & conf_, int a, triqs::gf::gf_view<triqs::gf::legendre> & Gl_):
  BaseType(), conf(conf_), a_level(a), Gl(Gl_), name(to_string("G(legendre)",a)){Gl()=0;}

 void accumulate(std::complex<double> signe);

 void collect_results( boost::mpi::communicator const & c){
   BaseType::collect_results(c);
   double Z_qmc ( real(this->acc_sign));

   for (auto l : Gl.mesh()) {
     Gl(l) = -(sqrt(2.0*l+1.0)/(Z_qmc*conf.Beta)) * Gl(l);
   }

   auto res = triqs::make_clone(Gl);
   auto g_loc = triqs::make_clone(Gl);
   boost::mpi::reduce(c, g_loc, res, std::plus<triqs::gf::gf<triqs::gf::legendre>>(),0);
   boost::mpi::broadcast(c,res,0);
   Gl = res;

 }

};

#endif
