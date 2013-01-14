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
#ifndef TRIQS_GF_LOCAL_FOURIER_MATSU_H 
#define TRIQS_GF_LOCAL_FOURIER_MATSU_H

#include <triqs/gf/imfreq.hpp> 
#include <triqs/gf/imtime.hpp> 

namespace triqs { namespace gf { 

 // First the implementation of the fourier transform
 void fourier_impl         (gf_view<imfreq> &gw , gf_view<imtime> const & gt);
 void inverse_fourier_impl (gf_view<imtime> &gt,  gf_view<imfreq> const & gw);

 // Then a good old function make a new gf 
 /*
 gf<imfreq> fourier (gf_view<imtime> const & gt) { 
  auto gw = imfreq::make_gf(gt.domain().beta, gt.domain().statistic,gt.data_view().shape().pop(),gt.mesh().size(), gt(freq_infty()));
  auto V = gw();
  fourier_impl(V,gt);
  return gw;
 }

 gf<imtime> inverse_fourier (gf_view<imfreq> const & gw) { 
  auto gt = imtime::make_gf(gw.domain().beta, gw.domain().statistic,gw.data_view().shape().pop(),gw.mesh().size(), gw(freq_infty()));
  auto V = gt();
  inverse_fourier_impl(V,gw);
  return gt;
 }
 */

 // Finally the lazy system for the = operator for views....
 namespace tags { struct fourier{}; }

 gf_keeper<tags::fourier,imtime> lazy_fourier         (gf_view<imtime> const & g) { return g;}
 gf_keeper<tags::fourier,imfreq> lazy_inverse_fourier (gf_view<imfreq> const & g) { return g;}

 void triqs_gf_view_assign_delegation( gf_view<imfreq> &g, gf_keeper<tags::fourier,imtime> const & L) { fourier_impl (g,L.g);}
 void triqs_gf_view_assign_delegation( gf_view<imtime> &g, gf_keeper<tags::fourier,imfreq> const & L) { inverse_fourier_impl(g,L.g);}

}}
#endif

