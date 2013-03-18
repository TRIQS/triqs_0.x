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
#include <fftw3.h>

namespace triqs { namespace gf { 

  namespace impl_local_real {
    dcomplex I(0,1);
    double pi = std::acos(-1);
    inline dcomplex th_expo(double t, double a ) { return (t < 0 ? 0 : -I * exp(-a*t)); }
    inline dcomplex th_expo_neg(double t, double a ) { return (t > 0 ? 0 : I * exp(a*t)); }
  }

  //--------------------------------------------------------------------------------------

  void fourier_impl(gf_view<refreq> & gw, gf_view<retime> const & gt) {

    using namespace impl_local_real;

    size_t L;
    switch(gt.mesh().kind()) {
      case half_bins: L = gt.mesh().size(); break;
      case full_bins: L = gt.mesh().size()-1; break;
      case without_last: L = gt.mesh().size(); break;
    }

    if (gw.mesh().size() != gt.mesh().size()) TRIQS_RUNTIME_ERROR << "Meshes are different";
    double test = std::abs(gt.mesh().delta() * gw.mesh().delta() * L / (2*pi) -1);
    if (test > 1.e-10) TRIQS_RUNTIME_ERROR << "Meshes are not compatible";

    const double tmin = gt.mesh().x_min() + (gt.mesh().kind() == half_bins ? 0.5 : 0.0) * gt.mesh().delta();
    const double wmin = gw.mesh().x_min() + (gw.mesh().kind() == half_bins ? 0.5 : 0.0) * gw.mesh().delta();

    auto ta = gt(freq_infty());
    tqa::vector<dcomplex> g_in(gt.mesh().size()), g_out(gw.mesh().size());

    for (size_t n1=0; n1<gw.data_view().shape()[0];n1++) {
      for (size_t n2=0; n2<gw.data_view().shape()[1];n2++) {

        dcomplex t1 = ta(1)(n1,n2), t2= ta(2)(n1,n2);
        dcomplex a1 = (t1 - I * t2)/2, a2 = (t1 + I * t2)/2;

        g_in() = 0;

        for (auto & t : gt.mesh()) {
          g_in(t.index()) = (gt(t)(n1,n2) - (a1*th_expo(t,1) + a2*th_expo_neg(t,1))) * std::exp(I*t*wmin);
        }

        details::fourier_base(g_in, g_out, L, true);

        for (auto & w : gw.mesh()) {
          gw(w)(n1,n2) = gt.mesh().delta() * std::exp(I*w*tmin) * std::exp(-I*wmin*tmin) * g_out(w.index()) + (a1/(w+I) + a2/(w-I));
        }
      }
    }

    // set tail
    gw.singularity_view() = gt.singularity_view();

  }

  //---------------------------------------------------------------------------

  void inverse_fourier_impl (gf_view<retime> & gt, gf_view<refreq> const & gw) { 

    using namespace impl_local_real;

    size_t L;
    switch(gw.mesh().kind()) {
      case half_bins: L = gw.mesh().size(); break;
      case full_bins: L = gw.mesh().size()-1; break;
      case without_last: L = gw.mesh().size(); break;
    }

    if (gw.mesh().size() != gt.mesh().size()) TRIQS_RUNTIME_ERROR << "Meshes are different";
    double test = std::abs(gt.mesh().delta() * gw.mesh().delta() * L / (2*pi) -1);
    if (test > 1.e-10) TRIQS_RUNTIME_ERROR << "Meshes are not compatible";

    const double tmin = gt.mesh().x_min() + (gt.mesh().kind() == half_bins ? 0.5 : 0.0) * gt.mesh().delta();
    const double wmin = gw.mesh().x_min() + (gw.mesh().kind() == half_bins ? 0.5 : 0.0) * gw.mesh().delta();

    auto ta = gw(freq_infty());
    tqa::vector<dcomplex> g_in(gw.mesh().size()), g_out(gt.mesh().size());

    for (size_t n1=0; n1<gt.data_view().shape()[0];n1++) {
      for (size_t n2=0; n2<gt.data_view().shape()[1];n2++) {

        dcomplex t1 = ta(1)(n1,n2), t2 = ta(2)(n1,n2);
        dcomplex a1 = (t1 - I * t2)/2, a2 = (t1 + I * t2)/2;
        g_in() = 0;

        for (auto & w: gw.mesh()) {
          g_in(w.index()) = (gw(w)(n1,n2) - (a1/(w+I) + a2/(w-I))) * std::exp(-I*w*tmin);
        }

        details::fourier_base(g_in, g_out, L, false);

        const double corr = 1.0/(gt.mesh().delta()*L);
        for (auto & t : gt.mesh()) {
          gt(t)(n1,n2) = corr * std::exp(I*wmin*tmin) * std::exp(-I*wmin*t) *
                                g_out(t.index()) + (a1*th_expo(t,1) + a2*th_expo_neg(t,1));
        }
      }
    }

    // set tail
    gt.singularity_view() = gw.singularity_view();

  }

  //---------------------------------------------------------------------------

  gf_keeper<tags::fourier,retime> lazy_fourier         (gf_view<retime> const & g) { return g;}
  gf_keeper<tags::fourier,refreq> lazy_inverse_fourier (gf_view<refreq> const & g) { return g;}

  void triqs_gf_view_assign_delegation( gf_view<refreq> &g, gf_keeper<tags::fourier,retime> const & L) { fourier_impl (g,L.g);}
  void triqs_gf_view_assign_delegation( gf_view<retime> &g, gf_keeper<tags::fourier,refreq> const & L) { inverse_fourier_impl(g,L.g);}

}}
