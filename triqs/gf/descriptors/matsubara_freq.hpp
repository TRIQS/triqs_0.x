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
#ifndef TRIQS_GF_MATSUBARA_FREQ_H
#define TRIQS_GF_MATSUBARA_FREQ_H
#include "../tools.hpp"
#include "../gf.hpp"
#include "../local/tail.hpp"
#include "../gf_proto.hpp"

namespace triqs { namespace gf { 

 struct matsubara_freq { 

  /// A tag to recognize the function 
  struct tag {};

  /// The domain
  struct domain_t {
   typedef std::complex<double> point_t;
   double beta;
   statistic_enum statistic;
   domain_t (double Beta, statistic_enum s = Fermion): beta(Beta), statistic(s){}
   bool operator == (domain_t const & D) { return ((std::abs(beta - D.beta)<1.e-15) && (statistic == D.statistic));}
  };

  /// The Mesh
  struct mesh_t {
   typedef matsubara_freq::domain_t domain_t;
   typedef size_t index_t;

   mesh_t (double Beta=1, statistic_enum s=Fermion, size_t N_max=1025): _dom(Beta,s), n_max_(N_max), pi_over_beta(std::acos(-1)/Beta), sh(s==Fermion? 1:0) {}

   domain_t const & domain() const { return _dom;}

   size_t size() const{ return n_max_;}

   // Conversions point <-> index <-> linear_index
   std::complex<double> index_to_point(index_t n) const {return std::complex<double>(0,(2*n+sh)*pi_over_beta);}
   size_t index_to_linear(index_t n) const { return n;}
   
   /// The wrapper for the mesh point
   struct mesh_point_t : arith_ops_by_cast<mesh_point_t,  std::complex<double> > {
    mesh_t const & m;  index_t index; mesh_point_t( mesh_t const & mesh, index_t const & index_=0): m(mesh), index(index_) {}
    void advance() { ++index;}
    typedef std::complex<double> cast_t;
    operator cast_t() const { return m.index_to_point(index);} 
   };

   // Accessing a point of the mesh
   mesh_point_t operator[](index_t i) const { return mesh_point_t(*this,i);}
   
   // Iterating on all the points...
   typedef  mesh_pt_generator<mesh_t> iterator;
   iterator begin() const { return iterator (*this);}
   iterator end()   const { return iterator (*this, true);}

   /// Mesh comparison
   friend bool operator == (mesh_t M1, mesh_t M2) { return ((M1._dom == M2._dom) && (M1.n_max_ ==M2.n_max_) );}

   private:
   domain_t _dom;
   size_t n_max_; 
   double pi_over_beta;
   index_t sh;
  };//end mesh_t

  /// The Auxiliary data
  typedef local::tail   auxiliary_data_t;

  /// The target
  typedef arrays::matrix<std::complex<double> >     target_t;
  //typedef arrays::matrix<std::complex<double>, arrays::Option::Fortran >     target_t;
  typedef typename target_t::view_type            target_view_t;

  /// Arity (number of argument in calling the function)
  static const int arity =1;

  /// All the possible calls of the gf
  template<typename D, typename T>
   target_view_t operator() (mesh_t const & mesh, D const & data, T const & t, long  n)  const {return data(arrays::range(), arrays::range(),n); } 

  template<typename D, typename T>
   local::tail_view operator()(mesh_t const & mesh, D const & data, T const & t, freq_infty const &) const {return t;} 

  /// How to fill a gf from an expression (RHS)
  template<typename D, typename T, typename RHS> 
   static void assign_from_expression (mesh_t const & mesh, D & data, T & t, RHS rhs) { 
    // access to the data . Beware, we view it as a *matrix* NOT an array... (crucial for assignment to scalars !) 
    for (size_t u=0; u<mesh.size(); ++u)  { target_view_t( data(tqa::range(),tqa::range(),u)) = rhs(mesh[u]); }
    t = rhs( local::tail::omega(t.shape(),t.size()));
    // if f is an expression, replace the placeholder with a simple tail. If f is a function callable on freq_infty, 
    // it uses the fact that tail_non_view_t can be casted into freq_infty 
   }

  /// Name of the auxiliary data e.g. in hdf5
  static std::string auxiliary_data_name() { return "tail";}

 };

 // -------------------------------   Expression template --------------------------------------------------

 // A trait to identify objects that have the concept ImmutableGfMatsubaraFreq
 template<typename G> struct ImmutableGfMatsubaraFreq : boost::is_base_of<typename matsubara_freq::tag,G> {};  

 // This defines the expression template with boost::proto (cf gf_proto.hpp).
 TRIQS_GF_DEFINE_OPERATORS(matsubara_freq,local::is_scalar_or_element,ImmutableGfMatsubaraFreq);

}}
#endif

