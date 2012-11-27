
/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012 by M. Ferrero, O. Parcollet, I. Krivenko
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


// TO DO : move to new GF. PUT in NAMESAPCE . clean array

#include "pade.hpp"
#include <gmpxx.h>
//#include<triqs/arrays.hpp>

typedef std::complex<double> dcomplex;

// This implementation is based on a Fortran code written by
// A. Poteryaev <Alexander.Poteryaev _at_ cpht.polytechnique.fr>
//
// The original algorithm is described in
// H. J. Vidberg, J. W. Serene. J. Low Temp. Phys. 29, 3-4, 179 (1977)

struct gmp_complex {
 mpf_class re, im;
 gmp_complex operator* (const gmp_complex &rhs){ return { rhs.re*re-rhs.im*im, rhs.re*im+rhs.im*re }; }
 friend gmp_complex inverse (const gmp_complex &rhs){ mpf_class d=rhs.re*rhs.re + rhs.im*rhs.im; return { rhs.re/d, -rhs.im/d}; }
 gmp_complex operator/ (const gmp_complex &rhs){ return (*this)*inverse(rhs); }
 gmp_complex operator+ (const gmp_complex &rhs){ return { rhs.re + re, rhs.im + im }; }
 gmp_complex operator- (const gmp_complex &rhs){ return { re - rhs.re, im - rhs.im }; }
 friend mpf_class real(const gmp_complex &rhs) { return rhs.re; }
 friend mpf_class imag(const gmp_complex &rhs) { return rhs.im; }
 gmp_complex& operator= (const std::complex<double> &rhs) { re = real(rhs); im = imag(rhs); return *this; }
 friend std::ostream & operator << (std::ostream & out,gmp_complex const & r) { return out << " gmp_complex("<<r.re<<","<<r.im<<")"<<std::endl ;}
};


class Pade_approximant {

 Array<dcomplex,1> z_in; // Input complex frequency points
 Array<dcomplex,1> a; // Pade coefficients

 //arrays::vector<dcomplex> z_in; // Input complex frequency points
 //arrays::vector<dcomplex> a; // Pade coefficients

 public:

 static const int GMP_default_prec = 256;    // Precision of GMP floats to use during a Pade coefficients calculation.
 Pade_approximant(const Array<dcomplex,1> & z_in_, const Array<dcomplex,1> & u_in) : z_in(z_in_) {
  //Pade_approximant(const arrays::vector<dcomplex> & z_in_, const arrays::vector<dcomplex> & u_in) : z_in(z_in_) {
  int N = z_in.size();
  a.resize(N);

  // Change the default precision of GMP floats.
  // MF: I have put an unsigned long instead of mp_bitcnt_t here so that it's compatible with GMP >= 4.3
  unsigned long old_prec = mpf_get_default_prec();
  mpf_set_default_prec(GMP_default_prec);  // How do we determine it?

  //arrays::array<gmp_complex,2> g(N,N);
  Array<gmp_complex,2> g(N,N);
  gmp_complex MP_0 = {0.0, 0.0};
  g = MP_0;
  for (int f = 0; f<N; ++f) { g(0,f) = u_in(f); };

  gmp_complex MP_1 = {1.0, 0.0};

  for(int p=1; p<N; ++p)
   for(int j=p; j<N; ++j){
    gmp_complex x = g(p-1,p-1)/g(p-1,j) - MP_1;
    gmp_complex y; y = z_in(j)-z_in(p-1);
    g(p,j) = x/y;
   }

  for(int j=0; j<N; ++j){
   gmp_complex gj(g(j,j));
   a(j) = dcomplex(real(gj).get_d(),imag(gj).get_d());
  }

  // Restore the precision.
  mpf_set_default_prec(old_prec);
 }

 dcomplex operator()(dcomplex e) const
 {
  dcomplex A1(0);
  dcomplex A2 = a(0);
  dcomplex B1(1.0), B2(1.0);

  int N = a.size();
  for(int i=0; i<=N-2; ++i){
   dcomplex Anew = A2 + (e - z_in(i))*a(i+1)*A1;
   dcomplex Bnew = B2 + (e - z_in(i))*a(i+1)*B1;
   A1 = A2; A2 = Anew;
   B1 = B2; B2 = Bnew;
  }

  return A2/B2;
 }

 };


 void pade(GF_Bloc_ImFreq const & Gw, GF_Bloc_ReFreq & Ge, int N_Matsubara_Frequencies, double Freq_Offset)
 {
  check_have_same_structure (Gw,Ge,false,true);
  assert (Gw.mesh.index_min==0);
  assert (Ge.mesh.index_min==0);

  double Beta = Gw.Beta;
  double omegaShift = (Gw.Statistic==Fermion ? 1 : 0);

  Array<dcomplex,1> z_in(N_Matsubara_Frequencies);
  firstIndex i;
  z_in = I*Pi/Beta*(2*i+omegaShift);

  // Just copy the tail. It doesn't need to conform to the Pade approximant.
  Gw.tail = Ge.tail;

  int N1 = Gw.N1, N2 = Gw.N2;
  for (int n1=1; n1<=N1;n1++)
   for (int n2=1; n2<=N2;n2++){
    int N = Gw.mesh.len();

    Array<dcomplex,1> u_in(N_Matsubara_Frequencies);     // Values at the Matsubara frequencies
    Array<dcomplex,1> a(N_Matsubara_Frequencies);        // Pade coefficients

    for(int i=0; i < N_Matsubara_Frequencies; ++i){
     u_in(i) = (i < N ? Gw.data_const(n1,n2,i) : Gw.tail.eval(z_in(i))(n1,n2));
    }

    Pade_approximant PA(z_in,u_in);

    int Ne = Ge.mesh.len();
    Ge.zero();
    for (int i=0; i < Ne; ++i) {
     dcomplex e = Ge.mesh[i] + I*Freq_Offset;
     Ge.data(n1,n2,i) = PA(e);
    }
   }
 }
