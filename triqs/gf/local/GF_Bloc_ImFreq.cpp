/
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

#include "GF_Bloc_ImFreq.hpp"

const bool Green_Function_Are_Complex_in_time = false;

GF_Bloc_ImFreq::GF_Bloc_ImFreq (boost::python::object IndicesL_,
				boost::python::object IndicesR_,
				PyObject * Data,
				boost::shared_ptr<MeshGF> Mesh,
				boost::shared_ptr<TailGF> Tail):
  GF_Bloc_Base<COMPLEX>(IndicesL_,IndicesR_, Data,Mesh,Tail)
{}
  
//------------------------------------------------------------------------------------------------

GF_Bloc_ImFreq::GF_Bloc_ImFreq(const GF_Bloc_ImFreq & Gin):
  GF_Bloc_Base<COMPLEX>(Gin)
{}

//------------------------------------------------------------------------------------------------

inline COMPLEX F(COMPLEX a,double b,double Beta) {return -a/(1+exp(-Beta*b));}

PyArray<COMPLEX,2> GF_Bloc_ImFreq::density() const
{
  if (!tail.is_decreasing_at_infinity())  TRIQS_RUNTIME_ERROR<<" GF_Bloc_ImFreq::density : Green Function is not as 1/omega or less !!!";
  
  Array<COMPLEX,2> dens_part(N1,N2,fortranArray);
  Array<COMPLEX,2> dens_tail(N1,N2,fortranArray);
  PyArray<COMPLEX,2> dens(N1,N2,FortranOrder); //dens(Nnd,Nnd,fortranArray);
  dens_part=0;dens=0;dens_tail=0;
  
  for (int n1=1; n1<=N1;n1++) 
    for (int n2=1; n2<=N2;n2++)
      {
	COMPLEX d=tail[1](n1,n2) , A=tail[2](n1,n2),B = tail[3](n1,n2) ;
	double b1 = 0,b2 =1, b3 =-1;
	COMPLEX a1 = d-B, a2 = (A+B)/2, a3 = (B-A)/2,om1;
	for (int om = mesh.index_min; om<=mesh.index_max; om++)
	  {
	    om1 = mesh[om];
	    dens_part(n1,n2)+= data(n1,n2,om)  -  (a1/(om1 - b1) + a2 / (om1-b2) + a3/(om1-b3)); 
	  }
	// If  the Green function are NOT complex, then one use the symmetry property
	// fold the sum and get a factor 2
	//double fact = (Green_Function_Are_Complex_in_time ? 1 : 2);
	//dens_part(n1,n2) = dens_part(n1,n2)*(fact/Beta)  + (d + F(a1,b1,Beta) + F(a2,b2,Beta)+ F(a3,b3,Beta));
	//if (!Green_Function_Are_Complex_in_time) dens_part  = 0+real(dens_part);
	dens_part(n1,n2) = dens_part(n1,n2)/Beta;
	dens_tail(n1,n2) = d + F(a1,b1,Beta) + F(a2,b2,Beta)+ F(a3,b3,Beta);
      }
   
  for (int n1=1; n1<=N1;n1++)
    for (int n2=n1; n2<=N2;n2++)
      {
  	dens_part(n1,n2) = dens_part(n1,n2) + real(dens_part(n2,n1)) - I * imag(dens_part(n2,n1)) + dens_tail(n1,n2);
  	dens_part(n2,n1) = real(dens_part(n1,n2)) - I * imag(dens_part(n1,n2));
      }
  

  (Array<COMPLEX,2>)dens = dens_part; 
  return(dens);
}

/*
  substraction : 
   a1/(om - b1) + a2/(om - b2) + a3/(om - b3) 

   expansion is : 
    1/om :   a1 + a2 + a3 = d
    1/om^2 : a1 b1 + a2 b2 + a3 b3 = A
    1/om^3 : a1 b1 + a2 b2^2 + a3 b3^2 = B

  Fourier inverse is : 
    \sum_i  a_i F_b(tau)
  with F_b(tau) = - exp(-b tau) /(1+ exp(-beta b))
  
  F_b(tau = 0^-) = F_b(0^+) + disc = disc - 1/(1 + exp( - beta* b)) 

*/
