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

#include <triqs/gf/matsubara_freq.hpp> 
#include <triqs/gf/matsubara_time.hpp> 

namespace triqs { namespace gf { 

 /*template<typename Descriptor, typename T,
 class fourier_lazy {
  typedef view_type_if_exists_else_type<T>::type keeper;
  mutable std::shared_ptr < gf<Descriptor> > _id;
  
  fourier_lazy( T const & x) : keeper(x){}

  gf<Descriptor> const & id() const { if (!_id) _id = std::make<gf<Descriptor>>(...);}

  typedef typename gf<Descriptor>::mesh_t             mesh_t;
  typedef typename gf<Descriptor>::singularity_view_t singularity_view_t;
  typedef typename gf<Descriptor>::data_view_t        data_view_t;

  mesh_t const & mesh() const                 { return id()->mesh();}
  const data_view_t data_view() const         { return id()->data_view();}
  const singularity_view_t singularity_view() const { return id()->singularity_view();}
  }; */


 template<typename Tag, typename D> struct gf_lazy_impl{ gf_view<D> g; gf_lazy_impl (gf_view<D> const & g_) : g(g_) {} };

 namespace tags { struct fourier{}; }

 void triqs_gf_view_assign_delegation( gf_view<matsubara_freq> &g, gf_lazy_impl<tags::fourier,matsubara_time> const & L) { fourier_direct_impl (g,L.g);}
 void triqs_gf_view_assign_delegation( gf_view<matsubara_time> &g, gf_lazy_impl<tags::fourier,matsubara_freq> const & L) { fourier_inverse_impl(g,L.g);}
 
 gf_lazy_impl<tags::fourier,matsubara_time> lazy_fourier_direct  (gf_view<matsubara_time> const & g) { return g;}
 gf_lazy_impl<tags::fourier,matsubara_freq> lazy_fourier_inverse (gf_view<matsubara_freq> const & g) { return g;}

 void fourier_direct_impl  (gf_view<matsubara_freq> &gw , gf_view<matsubara_time> const & gt);
 void fourier_inverse_impl (gf_view<matsubara_time> &gt,  gf_view<matsubara_freq> const & gw);

 gf<matsubara_freq> fourier_direct (gf_view<matsubara_time> const & gt) { 
  auto gw = matsubara_freq::make_gf(gt.domain().beta, gt.domain().statistic,gt.mesh().size(), gt.shape(),gt(freq_infty()));
  auto V = gw();
  fourier_direct_impl(V,gt);
  return gw;
 }

 gf<matsubara_time> fourier_inverse (gf_view<matsubara_freq> const & gw) { 
  auto gt = matsubara_time::make_gf(gw.domain().beta, gw.domain().statistic,gw.mesh().size(), gw.shape(),gw(freq_infty()));
  auto V = gt();
  fourier_inverse_impl(V,gw);
  return gt;
 }


}}
#endif


