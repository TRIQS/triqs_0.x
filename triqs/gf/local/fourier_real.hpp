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
#ifndef TRIQS_GF_LOCAL_FOURIER_REAL_H 
#define TRIQS_GF_LOCAL_FOURIER_REAL_H

#include "fourier_base.hpp"
#include <triqs/gf/refreq.hpp> 
#include <triqs/gf/retime.hpp> 

namespace triqs { namespace gf { 

 // First the implementation of the fourier transform
 void fourier_impl         (gf_view<refreq> &gw , gf_view<retime> const & gt);
 void inverse_fourier_impl (gf_view<retime> &gt,  gf_view<refreq> const & gw);

 // Then a good old function make a new gf 
 /*
 gf<refreq> fourier (gf_view<retime> const & gt) { 
  auto gw = refreq::make_gf(gt.domain().beta, gt.domain().statistic,gt.data_view().shape().pop(),gt.mesh().size(), gt(freq_infty()));
  auto V = gw();
  fourier_impl(V,gt);
  return gw;
 }

 gf<retime> inverse_fourier (gf_view<refreq> const & gw) { 
  auto gt = retime::make_gf(gw.domain().beta, gw.domain().statistic,gw.data_view().shape().pop(),gw.mesh().size(), gw(freq_infty()));
  auto V = gt();
  inverse_fourier_impl(V,gw);
  return gt;
 }
 */

 gf_keeper<tags::fourier,retime> lazy_fourier         (gf_view<retime> const & g);
 gf_keeper<tags::fourier,refreq> lazy_inverse_fourier (gf_view<refreq> const & g);

 void triqs_gf_view_assign_delegation( gf_view<refreq> &g, gf_keeper<tags::fourier,retime> const & L);
 void triqs_gf_view_assign_delegation( gf_view<retime> &g, gf_keeper<tags::fourier,refreq> const & L);

}}
#endif

