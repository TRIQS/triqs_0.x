/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012 by M. Ferrero, O. Parcollet
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

#ifndef TRIQS_GF_LOCAL_DENSITY_H
#define TRIQS_GF_LOCAL_DENSITY_H 

#include <triqs/lazy/core.hpp>
#include <triqs/arrays/array.hpp>
#include <triqs/arrays/matrix.hpp>

namespace triqs { namespace gf { namespace local {

 template< typename GfType, typename Enable =void > struct density_impl; 
 template< typename GfType> typename GfType::result_type density(GfType const & x) { return density_impl<GfType>::invoke(x);}

}}}

namespace triqs { namespace lazy {
 TRIQS_LAZY_MAKE_FNT_LAZY (1, density);
}}

// ------------- Implementations -----------------------

namespace triqs { namespace gf { namespace local {

 //-------------------------------------------------------
 // For Imaginary Matsubara Frequency functions
 // ------------------------------------------------------
 template< typename GfType> struct density_impl<GfType, typename boost::enable_if< typename has_gf_loc_concept<GfType, domains::matsubara_freq> >::type > { 

  dcomplex F(dcomplex a,double b,double Beta) {return -a/(1+exp(-Beta*b));}

  static typename GfType::result_type invoke(GfType const & G) {

   typedef std::complex<double> dcomplex;

   if (!G.tail.is_decreasing_at_infinity())  TRIQS_RUNTIME_ERROR<<" GF_Bloc_ImFreq::density : Green Function is not as 1/omega or less !!!";

   BOOST_AUTO( sh = G.result_shape());
   const size_t N1=sh[0], N2 = sh[1];
   triqs::arrays::array<dcomplex,2> dens_part(sh), dens_tail(sh), dens(sh);
   dens_part=0;dens=0;dens_tail=0;

   for (size_t n1=0; n1<N1;n1++) 
    for (size_t n2=0; n2<N2;n2++)
    {
     dcomplex d= G.tail(1)(n1,n2) , A=G.tail(2)(n1,n2),B = G.tail(3)(n1,n2) ;
     double b1 = 0,b2 =1, b3 =-1;
     dcomplex a1 = d-B, a2 = (A+B)/2, a3 = (B-A)/2,om1;
     for (size_t om = mesh.index_min; om<=mesh.index_max; om++)
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

   for (size_t n1=0; n1<N1;n1++) 
    for (size_t n2=0; n2<N2;n2++)
    {
     dens_part(n1,n2) = dens_part(n1,n2) + real(dens_part(n2,n1)) - I * imag(dens_part(n2,n1)) + dens_tail(n1,n2);
     dens_part(n2,n1) = real(dens_part(n1,n2)) - I * imag(dens_part(n1,n2));
    }

   return dens_part;

  }
 }; 


#endif




