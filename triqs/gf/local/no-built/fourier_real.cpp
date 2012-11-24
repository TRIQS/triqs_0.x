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
#include "fourier_real.hpp"
#include "fourier_base.hpp"
#include <fftw3.h>

namespace triqs { namespace gf { 

 template<typename GfElementType> GfElementType convert_green ( dcomplex const & x) { return x;}
 template<> double convert_green<double> ( dcomplex const & x) { return real(x);}

 inline dcomplex TH_EXPO(double t, double a ) { return (t < 0 ? 0 : dcomplex(0,-1) * exp(-a *t));}
 inline dcomplex TH_EXPO_NEG(double t, double a ) { return (t > 0 ? 0 : dcomplex(0,1)  * exp(a *t));}

 inline dcomplex oneFermion(dcomplex a,double b,double tau,double Beta) {
  return -a*( b >=0 ? exp(-b*tau)/(1+exp(-Beta*b)) : exp(b*(Beta-tau))/(1+exp(Beta*b)) );
 }

 inline dcomplex oneBoson(dcomplex a,double b,double tau,double Beta) {
  return a*( b >=0 ? exp(-b*tau)/(exp(-Beta*b)-1) : exp(b*(Beta-tau))/(1-exp(b*Beta)) );
 }

 /*
 inline void test_conformity( gf<one_time>  const & gt, gf<freq> const & gw) {
  long Nt = gt.mesh().size();
  long Nw = gw.mesh().size();
  const double Deltat = gt.domain().t_max - gt.domain().t_min;
  const double fact = Deltat/Nt;
  double Pi = std::acos(-1);
  if (Nw!=Nt) 
   TRIQS_RUNTIME_ERROR<<" the grid of the two functions are not adapted : different number of point: time "<<Nt<<" frequency "<<Nw;

  auto dw = gw.mesh()[1] - gw.mesh()[0];
  double pi = std::acos(-1);
  double test = abs(Deltat* (dw) /(2*pi) -1); 
  if ( test > 1.e-10)
   TRIQS_RUNTIME_ERROR<<" the grid of the two functions are not adapted : Delta_t * Delta_omega != 2 pi\n"
    << " TimeMax - TimeMin = "<< Deltat << "\n Delta_omega = "<< dw 
    << "\n (TimeMax - TimeMin) *Delta_omega/ (2 pi)  = "<< test;
 }
*/
 //--------------------------------------------------------------------------------------

 gf<freq> fourier_direct (gf<one_time> const & gt) { 
  //test_conformity(gt,gw);

  double pi = std::acos(-1);
  dcomplex I(0,1);
  long Nt = gt.mesh().size();
  const double omega0 = 2*pi/(gt.domain().t_max - gt.domain().t_min);
  const double tmin = gt.domain().t_min;
  const double Deltat = gt.domain().t_max - gt.domain().t_min;
  const double Wmax = pi / Deltat * Nt;
  const double fact = Deltat/Nt;

  auto ta = gt(freq_infty());
  gf<freq> gw ( freq::mesh_t (Wmax,Nt), gt.shape(),ta); 
  tqa::vector<dcomplex> g_in(gt.mesh().size()), g_out (gw.mesh().size()); 

  for (int n1=0; n1<gw.shape()[0];n1++)
   for (int n2=0; n2<gw.shape()[1];n2++) {
    dcomplex t1= ta(1)(n1,n2), t2= ta(2)(n1,n2); 
    dcomplex a1 = (t1 - I * t2)/2, a2 = (t1 + I * t2)/2;
    g_in() = 0;
    for (auto & t : gt.mesh())  
     g_in(t.index) = gt(t)(n1,n2) - (  a1*TH_EXPO(t,1) + a2 * TH_EXPO_NEG(t,1) );

    details::fourier_base(g_in, g_out, true);

    for (auto & w : gw.mesh()) {
     int p = (w.index + Nt/2)%Nt;
     gw(w)(n1,n2) = g_out(p)*exp(I*om*(tmin+ fact/2))*fact + ( a1/(w + I) + a2 / (w - I)) ;
    }
   }
  return gw;
 }
 //---------------------------------------------------------------------------

 gf<one_time> fourier_inverse (gf<freq> const & gw) { 

  double pi = std::acos(-1);
  dcomplex I(0,1);
  long Nt = gt.mesh().size();
  long Nw = gw.mesh().size();
  const double tmin = gt.domain().t_min;
  const double Deltat = gt.domain().t_max - gt.domain().t_min;
  const double fact = Deltat/Nt;
  tqa::vector<dcomplex> g_in(gt.mesh().size()), g_out (gw.mesh().size()); 

  auto ta = gw(freq_infty());
  gf<one_time> gt ( one_time::mesh_t (Wmax,Nw), gw.shape(), ta);

  tqa::vector<dcomplex> g_in(gt.mesh().size()), g_out (gw.mesh().size()); 

  for (int n1=0; n1<gw.shape()[0];n1++)
   for (int n2=0; n2<gw.shape()[1];n2++) {
    dcomplex t1= ta(1)(n1,n2), t2= ta(2)(n1,n2);
    dcomplex a1 = (t1 - I * t2)/2, a2 = (t1 + I * t2)/2;
    g_in() = 0;

    for (auto & w : gw.mesh()) {
     int p = (w.index + Nt/2)%Nt;
     g_in(p) = gw(w)(n1,n2)/(exp(I*om*(tmin + fact/2))*fact) - ( a1/(om + I) + a2 / (om - I)) ;
    }

    details::fourier_base(g_in, g_out, false);

    const double corr =1.0/Nt;
    for (auto & t : gt.mesh())  
     gt(t)(n1,n2) = corr*( g_out(t.index) + (  a1*TH_EXPO(t,1) + a2 * TH_EXPO_NEG(t,1) ) );
   }

  return gt;
 }
}} // namespace
