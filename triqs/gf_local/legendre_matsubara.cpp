
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

#include "GF_Bloc_ImFreq.hpp"
#include "GF_Bloc_ImTime.hpp"
#include "GF_Bloc_ImLegendre.hpp"
#include "fourier.hpp"
#include <triqs/utility/legendre.hpp>

using namespace triqs::utility;

void legendre_matsubara_direct (GF_Bloc_ImLegendre const & Gl, GF_Bloc_ImFreq & Gw) {

  check_have_same_structure (Gw,Gl,false,true);

  Gw.zero();
  Gw.tail = Gl.tail;

  // Use the transformation matrix
  for (int n1=1; n1<=Gw.N1;n1++)
    for (int n2=1; n2<=Gw.N2;n2++)
      for (int om = Gw.mesh.index_min; om <= Gw.mesh.index_max; om++)
        for (int l = Gl.mesh.index_min; l <= Gl.mesh.index_max; l++)
          Gw.data(n1,n2,om) += legendre_T(om,l) * Gl.data_const(n1,n2,l);

}

void legendre_matsubara_inverse (GF_Bloc_ImFreq const & Gw, GF_Bloc_ImLegendre & Gl) {

  check_have_same_structure (Gw,Gl,false,true);
  Gl.zero();

  // Construct a temporary imaginary-time Green's function Gt
  // I set Nt time bins. This is ugly, one day we must code the direct
  // transformation without going through imaginary time
  const int Nt = 50000;
  PyArray<double,1> M(Nt,COrder);
  for (int i=0; i<Nt; i++) M(i) = Gw.Beta/(2.0*Nt)+i*Gw.Beta/double(Nt);
  boost::shared_ptr<MeshGF> MGFptr(new MeshGF(Imaginary_Legendre,Fermion,Gw.Beta,M.as_PyObjectPtrBorrowedRef()));
  boost::shared_ptr<TailGF> TGFptr(new TailGF(Gw.tail.OrderMin(),Gw.tail.OrderMax(),Gw.IndicesL,Gw.IndicesR));
  GF_Bloc_ImTime Gt(Gw.IndicesL, Gw.IndicesR, Py_None, MGFptr, TGFptr);

  // We first transform to imaginary time because it's been coded with the knowledge of the tails
  fourier_inverse (Gw, Gt, true);
  legendre_matsubara_inverse (Gt, Gl);

}


void legendre_matsubara_direct (GF_Bloc_ImLegendre const & Gl, GF_Bloc_ImTime & Gt) {

  Gt.zero();
  Gt.tail = Gl.tail;
  legendre_generator L;

  for (int i=Gt.mesh.index_min; i<=Gt.mesh.index_max; i++) {
    L.reset( 2*Gt.mesh[i].real()/Gt.Beta-1 );
    for (int j = Gl.mesh.index_min; j <= Gl.mesh.index_max; j++)
      Gt.data(Range::all(),Range::all(),i)+=sqrt(2*j+1)/Gt.Beta*Gl.data_const(Range::all(),Range::all(),j)*L.next();
  }

}

void legendre_matsubara_inverse (GF_Bloc_ImTime const & Gt, GF_Bloc_ImLegendre & Gl) {

  Gl.zero();
  legendre_generator L;

  // Do the integral over imaginary time
  for (int i=Gt.mesh.index_min; i<=Gt.mesh.index_max; i++) {
    L.reset( 2*Gt.mesh[i].real()/Gt.Beta-1 );
    for (int l = Gl.mesh.index_min; l <= Gl.mesh.index_max; l++)
      Gl.data(Range::all(),Range::all(),l) += sqrt(2*l+1) * L.next() * Gt.data_const(Range::all(),Range::all(),i);
  }
  Gl *= Gt.Beta/Gt.numberTimeSlices;

  // Now get the tail
  Gl.determine_tail();

}
