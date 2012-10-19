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
#ifndef TRIQS_GF_TIME_H
#define TRIQS_GF_TIME_H
#include "../tools.hpp"
#include "../gf.hpp"
#include "../local/tail.hpp"
#include "../gf_proto.hpp"

namespace triqs { namespace gf { 

 struct one_time {

  /// A tag to recognize the function 
  struct tag {};

  /// The domain
  struct domain_t {
   typedef double point_t; 
   double t_max; 
   domain_t (double t_max_) : t_max(t_max_) {} 
   bool operator == (domain_t const & D) { return ((std::abs(t_max - D.t_max)<1.e-15) );}
  };

  /// The Mesh
  struct mesh_t {
    
   typedef one_time::domain_t         domain_t;
   typedef size_t index_t; // again ...

   mesh_t (double t_max, size_t nt_max) : _dom(t_max),L(nt_max){}
   mesh_t () : _dom(1),L(0){}
   domain_t const & domain() const { return _dom;}
   size_t size() const {return L;}

   /// Conversions point <-> index <-> linear_index
   double  index_to_point (index_t ind) const {return ind * domain().t_max/L;} 
   size_t  index_to_linear(index_t ind) const {return ind;}   
   index_t point_to_index (double t) const {index_t res = { (size_t)floor(t/domain().t_max*L)}; return res;}    

   /// The wrapper for the mesh point
   struct mesh_point_t : arith_ops_by_cast<mesh_point_t,  double > {
    mesh_t const & m; index_t index; mesh_point_t( mesh_t const & m_, index_t const & index_=0): m(m_), index(index_) {};
    void advance() { ++index;}
    typedef double cast_t;
    operator cast_t() const { return m.index_to_point(index);}
   };

   /// Accessing a point of the mesh
   mesh_point_t operator[](index_t i) const { return mesh_point_t (*this,i);}

   /// Iterating on all the points...
   typedef  mesh_pt_generator<mesh_t> iterator;
   iterator begin() const { return iterator (*this);}
   iterator end()   const { return iterator (*this, true);}

   /// Mesh comparison
   friend bool operator == (mesh_t M1, mesh_t M2) { return ((M1._dom == M2._dom) && (M1.L ==M2.L) );}

   private:
   domain_t _dom;
   size_t L; 
  }; //end mesh_t

  /// The Auxiliary data
  typedef local::tail auxiliary_data_t;

  /// The target
  typedef arrays::matrix<std::complex<double> >     target_t;
//  typedef arrays::matrix<std::complex<double>, arrays::Option::Fortran >     target_t;
  typedef typename target_t::view_type                                       target_view_t;

  /// Arity (number of argument in calling the function)
  static const int arity =1;

  /// All the possible calls of the gf
  template<typename D, typename T>
   target_view_t operator() (mesh_t const & mesh, D const & data, T const & t, double t0)  const {
    return data(arrays::range(), arrays::range(),mesh.index_to_linear(mesh.point_to_index (t0))); 
   } 

  template<typename D, typename T>
   local::tail_view operator()(mesh_t const & mesh, D const & data, T const & t, freq_infty const &) const {return t;} 

  /// How to fill a gf from an expression (RHS)
  template<typename D, typename T, typename RHS> 
   static void assign_from_expression (mesh_t const & mesh, D & data, T & t, RHS rhs) { 
    for (size_t u=0; u<mesh.size(); ++u)  { target_view_t( data(tqa::range(),tqa::range(),u)) = rhs(mesh[u]); }
    t = rhs( local::tail::omega(t.shape(),t.size()));
   }
   
  /// Name of the auxiliary data e.g. in hdf5
  static std::string auxiliary_data_name() { return "tail";}

 };

 // -------------------------------   Expression template --------------------------------------------------

 // A trait to identify objects that have the concept ImmutableGfOneTime
 template<typename G> struct ImmutableGfOneTime : boost::is_base_of<typename one_time::tag,G> {};  

 // This defines the expression template with boost::proto (cf gf_proto.hpp).
 // TRIQS_GF_DEFINE_OPERATORS(times,local::is_scalar_or_element,ImmutableGfOneTime);

}}
#endif

