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
#ifndef TRIQS_GF_MATSUBARA_TIME_H
#define TRIQS_GF_MATSUBARA_TIME_H
#include "../tools.hpp"
#include "../gf.hpp"
#include "../local/tail.hpp"
#include "../gf_proto.hpp"

namespace triqs { namespace gf { 

 struct matsubara_time { 

  /// A tag to recognize the function 
  struct tag {};

  /// The domain
  struct domain_t {
   typedef double point_t;
   double beta;
   statistic_enum statistic;
   domain_t (double Beta, statistic_enum s = Fermion): beta(Beta), statistic(s){}
   bool operator == (domain_t const & D) const { return ((std::abs(beta - D.beta)<1.e-15) && (statistic == D.statistic));}
  };

  /// The Mesh
  struct mesh_t {

   typedef matsubara_time::domain_t domain_t;
   typedef size_t index_t;

   mesh_t (double Beta=1, statistic_enum s=Fermion, size_t N_time_slices=1025, double delta=0):
    _dom(Beta,s), n_time_slices_(N_time_slices),_delta(delta), beta_over_l(Beta/N_time_slices){}

   domain_t const & domain() const { return _dom;}

   size_t size() const{ return n_time_slices_;}
   double delta() const { return _delta;}

   /// Conversions point <-> index <-> linear_index
   double index_to_point(index_t n) const {return _delta + n * beta_over_l;}
   size_t index_to_linear(index_t n) const { return n;}
   //index_t point_to_index (domain_t::point_t const & tau) const {index_t res =(size_t) floor( (tau - _delta) / beta_over_l ); return res;} 

   /// The wrapper for the mesh point
   typedef mesh_point_d1<mesh_t> mesh_point_t;

   /// Accessing a point of the mesh
   mesh_point_t operator[](index_t i) const { return mesh_point_t(*this,i);}

   /// Iterating on all the points...
   typedef  mesh_pt_generator<mesh_t> iterator;
   iterator begin() const { return iterator (this);}
   iterator end()   const { return iterator (this, true);}

   /// Mesh comparison
   friend bool operator == (mesh_t const& M1, mesh_t const &M2) { return ((M1._dom == M2._dom) && (M1.n_time_slices_ ==M2.n_time_slices_) );}

   /// Write into HDF5
   friend void h5_write (tqa::h5::group_or_file fg, std::string subgroup_name, mesh_t const & m) {
    tqa::h5::group_or_file gr =  fg.create_group(subgroup_name);
    h5_write(gr,"Beta",m.domain().beta);
    h5_write(gr,"Statistic",(m.domain().statistic==Fermion ? 'F' : 'B'));
    h5_write(gr,"size",m.size());
    h5_write(gr,"delta",m.delta());
   }

   /// Read from HDF5
   friend void h5_read  (tqa::h5::group_or_file fg, std::string subgroup_name, mesh_t & m){
    tqa::h5::group_or_file gr = fg.open_group(subgroup_name);
    size_t Nmax; double delta, beta; char statistic;
    h5_read(gr,"Beta",beta);
    h5_read(gr,"Statistic",statistic);
    h5_read(gr,"size",Nmax);
    h5_read(gr,"delta",delta);
    m = mesh_t(beta,(statistic=='F' ? Fermion : Boson),Nmax,delta);
   }

   private:
   domain_t _dom;
   size_t n_time_slices_;
   double _delta;
   double beta_over_l;
  }; // end mesh_t

  /// The target
  typedef arrays::matrix<std::complex<double> >     target_t;
  //  typedef arrays::matrix<std::complex<double>, arrays::Option::Fortran >     target_t;
  typedef typename target_t::view_type            target_view_t;

  /// The tail
  typedef local::tail   singularity_t;

  /// Symmetry
  typedef nothing symmetry_t;

  /// Arity (number of argument in calling the function)
  static const int arity =1;

  /// All the possible calls of the gf
  //ERROR : give a double and interpolate
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

  static std::string h5_name() { return "matsubara_time";}
 };

 // -------------------------------   Expression template --------------------------------------------------

 // A trait to identify objects that have the concept ImmutableGfMatsubaraFreq
 template<typename G> struct ImmutableGfMatsubaraTime : boost::is_base_of<typename matsubara_time::tag,G> {};  

 // This defines the expression template with boost::proto (cf gf_proto.hpp).
 TRIQS_GF_DEFINE_OPERATORS(matsubara_time,local::is_scalar_or_element,ImmutableGfMatsubaraTime);

}}
#endif

