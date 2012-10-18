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
   typedef two_times::domain_t         domain_t;
   typedef typename domain_t::point_t  point_t;
   //typedef std::array<size_t,2> index_t; // again ...
   typedef tqa::mini_vector<size_t,2> index_t; // again ...

   mesh_t (double t_max, size_t n_max) : _dom(t_max),L(n_max){}
   mesh_t () : _dom(1),L(0){}

   domain_t const & domain() const { return _dom;}

   size_t size() const {return L*L;}

   // Conversions point <-> index <-> linear_index
   double index_to_point (index_t const & ind, zero_t) const {return ind[0] * domain().t_max/L;} 
   double index_to_point (index_t const & ind, one_t)  const {return ind[1] * domain().t_max/L;} 

   size_t index_to_linear(index_t const & ind) const { return ind[0] + ind[1]*L;}

   point_t index_to_point(index_t const & ind) const {
    double s= domain().t_max/L; domain_t::point_t res = {ind[0]*s, ind[1]*s}; return res;
   }
   
   index_t point_to_index(double t0, double t1) const {
    index_t res = { (size_t)floor(t0/domain().t_max*L), (size_t)floor(t1/domain().t_max*L)}; return res; 
   }    
   
   index_t point_to_index(point_t const & p) const { return point_to_index(p[0],p[1]);}

   /// The wrapper for the mesh point
   struct mesh_point_t { 
    mesh_t const & m; index_t index;

    /// The wrapper for one component of the mesh point
    template<typename Order> struct component : arith_ops_by_cast<component<Order>, double>  { 
     mesh_point_t const & p; 
     component (mesh_point_t const & p_):p(p_){}
     operator double() const { return p.m.index_to_point(p.index,Order());}
    };
    component<zero_t> t0;
    component<one_t>  t1;

    mesh_point_t(mesh_t const & m_, index_t index_ ) : m(m_), index(index_), t0(*this), t1(*this) { }
    mesh_point_t(mesh_t const & m_) : m(m_), index{{0,0}}, t0(*this), t1(*this) { }

    void advance() { ++index[0]; if (index[0]==m.L) { index[0]=0; ++index[1];}}
    typedef domain_t::point_t cast_t;
    operator cast_t() const { return m.index_to_point(index);}
   };

   // Accessing a point of the mesh
   mesh_point_t operator[](index_t i) const { return mesh_point_t (*this,i);}
   mesh_point_t operator()(size_t i, size_t j) const { index_t k= {i,j}; return mesh_point_t (*this, k);}

   // Iterating on all the points...
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
  typedef nothing auxiliary_data_t;

  /// The target
  typedef arrays::matrix<std::complex<double> >     target_t;
//  typedef arrays::matrix<std::complex<double>, arrays::Option::Fortran >     target_t;
  typedef typename target_t::view_type                                       target_view_t;

  /// Arity (number of argument in calling the function)
  static const int arity =2;

  /// All the possible calls of the gf
  template<typename D, typename T>
   target_view_t operator() (mesh_t const & mesh, D const & data, T const & t, double t1, double t2)  const {
    return data(arrays::range(), arrays::range(),mesh.index_to_linear(mesh.point_to_index (t1,t2)));
   } 

  /// How to fill a gf from an expression (RHS)
  template<typename D, typename T, typename RHS> 
   static void assign_from_expression (mesh_t const & mesh, D & data, T & t, RHS rhs) { 
    // access to the data . Beware, we view it as a *matrix* NOT an array... (crucial for assignment to scalars !) 
    for (auto & p : mesh) {  
     //std::cout  << " rhs "<< rhs(p.t0,p.t1)<< std::endl;
     //std::cout  << " p "<< p.index[0] << "  "<< p.index[1]<< std::endl;
     target_view_t( data(tqa::range(),tqa::range(),p.m.index_to_linear(p.index))) = rhs(p.t0,p.t1); }
    //for (size_t u=0; u<mesh.size(); ++u)  { target_view_t( data(tqa::range(),tqa::range(),u)) = rhs(mesh[u]); }
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

