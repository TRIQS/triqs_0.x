
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

#ifndef TRIQS_CTHYB1_BINNER_H
#define TRIQS_CTHYB1_BINNER_H

template< class GFType>
class gf_binner {
 GFType & G;
 const double L_over_Beta;
 public : 
 gf_binner ( GFType & G_): G(G_), L_over_Beta(1.0/G.mesh().delta()) {}

 // Binning operation
 void operator() (int alpha1, double t1, int alpha2, double t2,  typename GFType::data_t::value_type val)  { 
  if (t1>=t2) 
   G.data()(int(floor((t1 - t2)*L_over_Beta)), alpha1, alpha2) += val;
  else
   G.data()(int(floor((t1 - t2 + G.domain().beta)*L_over_Beta)), alpha1, alpha2) -=val;
 }

}; 


template< class GFType>
class gf_grid_evaluator {
 GFType const & G;
 const double L_over_Beta;
 public : 
 gf_grid_evaluator ( GFType const  & G_): G(G_), L_over_Beta(1.0/G.mesh().delta()) {}

 // For speed, it is better to make no interpolation and use a very fine grid.
 // otherwise the function is not really inlined.... and it is evaluated a lot !
 // alpha starts at ZERO
 typename GFType::data_t::value_type operator() (int alpha1, double tau1, int alpha2, double tau2 ) const 
 { 
  return ( (tau1>=tau2) ? G.data()(int(floor((tau1 - tau2)*L_over_Beta)), alpha1, alpha2) : 
    - G.data()(int(floor((tau1 - tau2 + G.domain().beta)*L_over_Beta)), alpha1, alpha2 ) );
 }

}; 

#endif

