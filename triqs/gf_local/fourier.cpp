
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

#include "GF_Bloc_ImTime.hpp"
#include "GF_Bloc_ImFreq.hpp"
#include "GF_Bloc_ReTime.hpp"
#include "GF_Bloc_ReFreq.hpp"
#include <fftw3.h>

inline COMPLEX oneFermion(COMPLEX a,double b,double tau,double Beta) {
 return -a*( b >=0 ? exp(-b*tau)/(1+exp(-Beta*b)) : exp(b*(Beta-tau))/(1+exp(Beta*b)) );
}

inline COMPLEX oneBoson(COMPLEX a,double b,double tau,double Beta) {
 return a*( b >=0 ? exp(-b*tau)/(exp(-Beta*b)-1) : exp(b*(Beta-tau))/(1-exp(b*Beta)) );
}
//--------------------------------------------------------------------------------------

void fourier_base(const Array<COMPLEX,1> &in, Array<COMPLEX,1> &out, bool direct) {
 
  // !!!! L must always be the number of time bins !!!!
  const int L( (direct ? in.extent(0) : out.extent(0)) );
  //const int L(max(in.extent(0),out.extent(0)));  <-- bug

  fftw_complex *inFFT, *outFFT;
  fftw_plan p;
  inFFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * L);
  outFFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * L);

  const COMPLEX * restrict in_ptr = in.dataFirst();
  COMPLEX * restrict out_ptr = out.dataFirst();
  const int imax(min(L,in.extent(0)));
  for (int i =0; i<imax; ++i) { inFFT[i][0] = real(in_ptr[i]); inFFT[i][1] = imag(in_ptr[i]);}
  for (int i =imax; i<L; ++i) { inFFT[i][0] = 0; inFFT[i][1] = 0;}

  p = fftw_plan_dft_1d(L, inFFT, outFFT, (direct ? FFTW_BACKWARD : FFTW_FORWARD), FFTW_ESTIMATE);
  fftw_execute(p); 

  const int jmax(min(L,out.extent(0)));
  for (int j =0; j<jmax; ++j)  {out_ptr[j] = COMPLEX(outFFT[j][0] ,outFFT[j][1]);}

  fftw_destroy_plan(p); 
  fftw_free(inFFT); fftw_free(outFFT);
}

//--------------------------------------------------------------------------------------
//
void fourier_direct  (GF_Bloc_ImTime const & Gt, GF_Bloc_ImFreq & Gw, bool time_mesh_starts_at_half_bin) { 
 double Beta = Gw.Beta;
 //assert(Gw.Statistic==Fermion); // Boson not implemented.
 check_have_same_structure (Gw,Gt,false,true);
 assert (Gw.mesh.index_min==0);

 Gw.zero();
 Array<COMPLEX,1>  G_out(Range(Gw.mesh.index_min, Gw.mesh.index_max)), 
  G_in(Range(Gt.mesh.index_min, Gt.mesh.index_max));//G_in(Gt.mesh.len());

 Gw.tail = Gt.tail;
 int N1 = Gt.N1, N2 = Gt.N2;

 double t;

 for (int n1=1; n1<=N1;n1++)
  for(int n2=1; n2<=N2;n2++)
  {
   COMPLEX d= Gw.tail[1](n1,n2), A= Gw.tail[2](n1,n2),B = Gw.tail[3](n1,n2);
   double b1, b2, b3;
   COMPLEX a1, a2, a3;
   if (Gw.Statistic == Fermion){ 
      b1 = 0; b2 =1; b3 =-1;
      a1 = d-B; a2 = (A+B)/2; a3 = (B-A)/2;
   }
   else {
      b1 = -0.5; b2 =-1; b3 =1;
      a1=4*(d-B)/3; a2=B-(d+A)/2; a3=d/6+A/2+B/3;
      }
   G_in = 0;

   for (int i =Gt.mesh.index_min; i <= Gt.mesh.index_max ; i++)
   {
    t = real(Gt.mesh[i]);
    if  (Gw.Statistic == Fermion){
    G_in(i) = Gt.data_const(n1,n2,i) - ( oneFermion(a1,b1,t,Beta) + oneFermion(a2,b2,t,Beta)+ oneFermion(a3,b3,t,Beta) );
    G_in(i) *= exp(I*Pi*t/Beta);
    }
    else {
     G_in(i) = Gt.data_const(n1,n2,i) - ( oneBoson(a1,b1,t,Beta) + oneBoson(a2,b2,t,Beta)+ oneBoson(a3,b3,t,Beta) );
    }
   }
   G_in *= Beta/Gt.numberTimeSlices;

   fourier_base(G_in, G_out, true);

   //	const COMPLEX C(I*Pi/Gt.mesh.lenTotal());
   for (int om = Gw.mesh.index_min; om<= Gw.mesh.index_max ; om++)
   {
    COMPLEX om1 = Gw.mesh[om];
    if (time_mesh_starts_at_half_bin) {
      Gw.data(n1,n2,om) = G_out(om)*exp(I*om*Pi/Gt.numberTimeSlices) + a1/(om1 - b1) + a2 / (om1-b2) + a3/(om1-b3);
    } else {
      Gw.data(n1,n2,om) = G_out(om) + a1/(om1 - b1) + a2 / (om1-b2) + a3/(om1-b3); 
    }
     

   }
  }
}

//---------------------------------------------------------------------------

static bool Green_Function_Are_Complex_in_time = false;
#define conj_green(x)  (x)
#define convert_green(x)  real(x)
#define green_to_complex(x) real_as_C(x)
void fourier_inverse (GF_Bloc_ImFreq const & Gw, GF_Bloc_ImTime & Gt, bool time_mesh_starts_at_half_bin) {

 if (abs(Gw.Beta - Gt.Beta)> 1.e-10) TRIQS_RUNTIME_ERROR<<"The time and frequency Green function are not at the same Beta";
 double Beta = Gw.Beta;
 //assert(Gw.Statistic==Fermion); // Boson not implemented.
 check_have_same_structure (Gw,Gt,false,true);
 assert (Gw.mesh.index_min==0);

 Gt.zero();
 Array<COMPLEX,1>  G_out(Range(Gt.mesh.index_min, Gt.mesh.index_max));
 Array<COMPLEX,1>  G_in(Range(Gw.mesh.index_min, Gw.mesh.index_max));

 Gt.tail = Gw.tail;
 int N1 = Gw.N1, N2 = Gw.N2;


 for (int n1=1; n1<=N1;n1++) 
  for (int n2=1; n2<=N2;n2++) 
  {
   COMPLEX d= Gt.tail[1](n1,n2), A= Gt.tail[2](n1,n2),B = Gt.tail[3](n1,n2);
   //double b1 = 0,b2 =1, b3 =-1;
   //COMPLEX a1 = d-B, a2 = (A+B)/2, a3 = (B-A)/2;
   double b1, b2, b3;
   COMPLEX a1, a2, a3;
   if (Gw.Statistic == Fermion){ 
      b1 = 0; b2 =1; b3 =-1;
      a1 = d-B; a2 = (A+B)/2; a3 = (B-A)/2;
   }
   else {
      b1 = -0.5; b2 =-1; b3 =1;
      a1=4*(d-B)/3; a2=B-(d+A)/2; a3=d/6+A/2+B/3;
      }
   G_in = 0;
   for (int i = Gw.mesh.index_min; i <= Gw.mesh.index_max ; i++) 
   {
    COMPLEX om1 = Gw.mesh[i];
    G_in(i) =  Gw.data_const(n1,n2,i) - ( a1/(om1 - b1) + a2 / (om1-b2) + a3/(om1-b3) );
    // you need an extra factor if the time mesh starts at 0.5*Beta/L
    if (time_mesh_starts_at_half_bin) G_in(i) *=  exp(-i*I*Pi/Gt.numberTimeSlices);
   }
   // for bosons GF(w=0) is divided by 2 to avoid counting it twice
   if (Gw.Statistic == Boson && !Green_Function_Are_Complex_in_time ) G_in(Gw.mesh.index_min) *= 0.5; 

   G_out = 0;
   fourier_base(G_in, G_out, false);

   // If  the Green function are NOT complex, then one use the symmetry property
   // fold the sum and get a factor 2
   double fact = (Green_Function_Are_Complex_in_time ? 1 : 2),t;
   G_out = G_out*(fact/Beta);
   for (int i =Gt.mesh.index_min; i<= Gt.mesh.index_max ; i++) 
   { 
    t= real(Gt.mesh[i]);
    if(Gw.Statistic == Fermion){
      Gt.data(n1,n2,i) = convert_green(G_out(i)*exp(-I*Pi*t/Beta) + oneFermion(a1,b1,t,Beta) + oneFermion(a2,b2,t,Beta)+ oneFermion(a3,b3,t,Beta) );}
    else {
      Gt.data(n1,n2,i) = convert_green(G_out(i) + oneBoson(a1,b1,t,Beta) + oneBoson(a2,b2,t,Beta)+ oneBoson(a3,b3,t,Beta) );}
   }
  }
}

//=============================================================================

inline COMPLEX TH_EXPO(double t, double a ) { return (t < 0 ? 0 : -I * exp(-a *t));}
inline COMPLEX TH_EXPO_NEG(double t, double a ) { return (t > 0 ? 0 : +I * exp(a *t));}

void test_conformity( GF_Bloc_ReTime const & Gt, GF_Bloc_ReFreq const & Gw) {
   const int N(Gw.mesh.len());
  const int Nt(Gt.mesh.len());
  const double Deltat((Gt.TimeMax - Gt.TimeMin));
  const double fact(Deltat/N);
  if (N!=Nt) 
   TRIQS_RUNTIME_ERROR<<" the grid of the two functions are not adapted : different number of point: time "<<Nt<<" frequency "<<N;

  double test = abs(Deltat* (Gw.mesh[1] - Gw.mesh[0]) /(2*pi) -1); 
  if ( test > 1.e-10)
   TRIQS_RUNTIME_ERROR<<" the grid of the two functions are not adapted : Delta_t * Delta_omega != 2 pi\n"
    << " TimeMax - TimeMin = "<< Deltat 
    << "\n Delta_omega = "<< (Gw.mesh[1] - Gw.mesh[0]) 
    << "\n (TimeMax - TimeMin) *Delta_omega/ (2 pi)  = "<< test;

}

void fourier_direct(GF_Bloc_ReTime const & Gt, GF_Bloc_ReFreq & Gw)
{
  assert(Gw.Statistic==Fermion); // Boson not implemented.
  check_have_same_structure (Gw,Gt,false,true);
  assert (Gw.mesh.index_min==0);
  assert (Gt.mesh.index_min==0);
  
  const int N(Gw.mesh.len());
  const int Nt(Gt.mesh.len());
  const double omega0(2*pi/(Gt.TimeMax - Gt.TimeMin));
  const double tmin( Gt.TimeMin);
  const double Deltat((Gt.TimeMax - Gt.TimeMin));
  const double fact(Deltat/N);
  test_conformity(Gt,Gw);

  Gw.zero();

  Array<COMPLEX,1>  G_out(N),G_in(N);
  double t; COMPLEX om;

  // check the frequency grid
  for (int i = 0; i <= N/2 ; i++) {
   om = Gw.mesh[i];
   assert ( abs (om - (-N/2 + i)*omega0) < 1.e-5);
  }
  for (int i = N/2+1; i < N ; i++) {
   om = Gw.mesh[i];
   assert ( abs (om - (-N/2  + i)*omega0) < 1.e-5);
  }

  Gw.tail = Gt.tail;
  int N1 = Gw.N1, N2 = Gw.N2;

  for (int n1=1; n1<=N1;n1++)
   for (int n2=1; n2<=N2;n2++)
   {
    COMPLEX t1= Gw.tail[1](n1,n2), t2= Gw.tail[2](n1,n2); // t3 = tail[3](n1,n2);
    COMPLEX a1 = (t1 - I * t2)/2, a2 = (t1 + I * t2)/2;
    G_in = 0;
    for (int i =0; i<N ; i++) {
     t= real(Gt.mesh[i]);
     G_in(i) =  Gt.data_const(n1,n2,i) - (  a1*TH_EXPO(t,1) + a2 * TH_EXPO_NEG(t,1) );
    }

    fourier_base(G_in, G_out, true);

    for (int i = 0; i < N ; i++) {
     om = Gw.mesh[i];
     int p = (i+ N/2)%N;
     Gw.data(n1,n2,i) = G_out(p)*exp(I*om*(tmin+ fact/2))*fact + ( a1/(om + I) + a2 / (om - I)) ;
    }

   }
}

//---------------------------------------------------------------------------

void fourier_inverse (GF_Bloc_ReFreq const & Gw, GF_Bloc_ReTime & Gt ) {
 assert(Gw.Statistic==Fermion); // Boson not implemented.
 check_have_same_structure (Gw,Gt,false,true);
 assert (Gw.mesh.index_min==0);

 const int N(Gw.mesh.index_max - Gw.mesh.index_min+1);
 const int Nt(Gt.mesh.index_max - Gt.mesh.index_min+1);
 const double tmin( Gt.TimeMin);
 const double Deltat((Gt.TimeMax - Gt.TimeMin));
 const double fact(Deltat/N);
 const double pi = acos(-1);
 test_conformity(Gt,Gw);

 Gt.zero();

 Array<COMPLEX,1>  G_out(N),G_in(N);
 double t; COMPLEX om;

 Gt.tail = Gw.tail;
 int N1 = Gt.N1, N2 = Gt.N2;

 for (int n1=1; n1<=N1;n1++) 
  for (int n2=1; n2<=N2;n2++) 
  {
   COMPLEX t1= Gt.tail[1](n1,n2), t2= Gt.tail[2](n1,n2); // t3 = tail[3](n1,n2);
   COMPLEX a1 = (t1 - I * t2)/2, a2 = (t1 + I * t2)/2;
   G_in = 0;
   for (int i = 0; i < N ; i++) {
    om = Gw.mesh[i];
    int p = (i+ N/2)%N;
    G_in(p) =  Gw.data_const(n1,n2,i)/(exp(I*om*(tmin + fact/2))*fact) - ( a1/(om + I) + a2 / (om - I)) ;
   }

   fourier_base(G_in, G_out, false);

   const double corr(1.0/N);
   for (int i =0; i<N ; i++) {
    t= real(Gt.mesh[i]);
    Gt.data(n1,n2,i)  = corr*( G_out(i) + (  a1*TH_EXPO(t,1) + a2 * TH_EXPO_NEG(t,1) ) );
   }
  }
}


