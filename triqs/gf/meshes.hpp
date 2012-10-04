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
#ifndef TRIQS_GF_LOCAL_MESHES_H
#define TRIQS_GF_LOCAL_MESHES_H

#include "domains.hpp"
#include "mesh_pt.hpp"

namespace triqs { namespace gf { namespace meshes { 

 //--------------------------------------------------------
 struct real_freq {};

 //--------------------------------------------------------
/*
 struct BZ_mesh{
  typedef domains::BZ domain_type;
  typedef long index_type;
  typedef arrays::range slice_args_type;

  domain_type const & domain() const {return dom_;}
  mesh_pt<BZ_mesh> operator[](index_type n) const  { return make_mesh_pt(*this,n);}
  domain_type::point_type embed(index_type const & n) const {return 0;}

  template<typename F> typename F::mv_type interpolate( F const & f, domain_type::point_type const & k) const {   return f((*this)[ 0 ]); } //just to fill in the fct

  protected:
  domain_type dom_;
 };
*/
}
}}

#endif




