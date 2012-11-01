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
#ifndef TRIQS_GF_TOOLS_H
#define TRIQS_GF_TOOLS_H
#include <triqs/utility/first_include.hpp>
#include <utility>
#include <boost/iterator/iterator_facade.hpp>
#include <triqs/utility/cint.hpp>
#include <triqs/clef/core.hpp>
#include <triqs/arrays/array.hpp>
#include <triqs/arrays/matrix.hpp>
#include <triqs/utility/proto/tools.hpp>
#include <triqs/arrays/h5/simple_read_write.hpp>
#include <triqs/arrays/proto/matrix_algebra.hpp>
#include <triqs/arrays/proto/array_algebra.hpp>
#include "triqs/utility/complex_ops.hpp"
#include <triqs/utility/view_tools.hpp>

namespace triqs { namespace gf {

 namespace tqa= triqs::arrays;
 namespace mpl=boost::mpl;
 namespace proto=boost::proto;
 namespace bf = boost::fusion; 
 namespace tup = triqs::utility::proto;

 //------------------------------------------------------

 typedef std::complex<double> dcomplex; 

 enum statistic_enum {Boson,Fermion};

 struct freq_infty{}; // the point at infinity

 //------------------------------------------------------

 struct nothing {
  template<typename... Args> nothing(Args...) {} // takes anything, do nothing..
  nothing() {}
  typedef void has_view_type_tag;     // Idiom : ValueView  
  typedef nothing view_type;
  typedef nothing non_view_type;
  void rebind (nothing){}
  void operator=(nothing) {}
  friend void h5_write (tqa::h5::group_or_file fg, std::string subgroup_name, nothing ) {}
  friend void h5_read  (tqa::h5::group_or_file fg, std::string subgroup_name, nothing ) {}
 }; 

 //------------------------------------------------------

}}
#endif
