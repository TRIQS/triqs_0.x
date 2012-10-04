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
#ifndef TRIQS_GF_2TIMES_H
#define TRIQS_GF_2TIMES_H
#include "../tools.hpp"
#include "../gf.hpp"
#include "../gf_proto.hpp"
#include <array>

namespace triqs { namespace gf { 

 struct two_times { 

  /// A tag to recognize the function 
  struct tag {};

  /// The domain
  struct domain_t {
   typedef std::array<double,2> point_t; // a fixed size (2) array of double... name is unfortunate, everything is called array !
   double t_max; 
   domain_t (double t_max_) : t_max(t_max_) {} 
   bool operator == (domain_t const & D) { return ((std::abs(t_max - D.t_max)<1.e-15) );}
  };

  /// The Mesh
  struct mesh_t {
   typedef two_times::domain_t domain_t;
   typedef std::array<size_t,2> index_t; // again ...

   mesh_t (double t_max, size_t n_max) : _dom(t_max),n_max_(n_max){}

   domain_t const & domain() const { return _dom;}

   size_t size() const{ return n_max_;}

   typedef mesh_pt<mesh_t, std::complex<double> > mesh_point_t;
   mesh_point_t operator[](index_t i) const { return mesh_point_t(*this,i);}

   std::complex<double> index_to_point(index_t n) const {return std::complex<double>(0,(2*n+sh)*pi_over_beta);}

   typedef  mesh_pt_generator<mesh_t> iterator;
   iterator begin() const { return iterator (*this);}
   iterator end()   const { return iterator (*this, true);}

   friend bool operator == (mesh_t M1, mesh_t M2) { return ((M1._dom == M2._dom) && (M1.n_max_ ==M2.n_max_) );}

   private:
   domain_t _dom;
   size_t n_max_; 
  };

  /// The Auxiliary data
  typedef nothing auxiliary_data_t;

  /// The target
  typedef arrays::matrix<std::complex<double>, arrays::Option::Fortran >     target_t;
  typedef typename target_t::view_type            target_view_t;

  /// Arity (number of argument in calling the function)
  static const int arity =2;

  /// All the possible calls of the gf
  template<typename D, typename T>
   target_view_t operator() (mesh_t const & mesh, D const & data, T const & t, double t1, double t2)  const {
    // compute index ?
    return data(arrays::range(), arrays::range(),n); 
   } 

  /// How to fill a gf from an expression (RHS)
  template<typename D, typename T, typename RHS> 
   static void assign_from_expression (mesh_t const & mesh, D & data, T & t, RHS rhs) { 
    // access to the data . Beware, we view it as a *matrix* NOT an array... (crucial for assignment to scalars !) 
    for (size_t u=0; u<mesh.size(); ++u)  { target_view_t( data(tqa::range(),tqa::range(),u)) = rhs(mesh[u]); }
   }

  /// Name of the auxiliary data e.g. in hdf5
  static std::string auxiliary_data_name() { return "tail";}

 };

 // -------------------------------   Expression template --------------------------------------------------

 // A trait to identify objects that have the concept ImmutableGfTwoTimes
 template<typename G> struct ImmutableGfTwoTimes : boost::is_base_of<typename two_times::tag,G> {};  

 // This defines the expression template with boost::proto (cf gf_proto.hpp).
 // TRIQS_GF_DEFINE_OPERATORS(two_times,local::is_scalar_or_element,ImmutableGfTwoTimes);

}}
#endif

