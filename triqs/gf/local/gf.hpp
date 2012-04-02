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
#ifndef TRIQS_GF_LOCAL_GF_H
#define TRIQS_GF_LOCAL_GF_H 

#include <triqs/lazy/core.hpp>
#include <triqs/arrays/array.hpp>
#include <triqs/arrays/matrix.hpp>
#include "./domains.hpp"

namespace triqs { namespace gf { namespace local {

 namespace tqa= triqs::arrays; namespace tql= triqs::lazy; namespace mpl= boost::mpl;
 using tqa::range;

 template<typename MeshType, bool IsView> class gf;

 namespace impl {

  template<typename LHS, typename RHS, typename Enable=void > struct assignment;

  struct nothing { // used in place of tail when there is none... an object that does .. nothing 
   typedef nothing view_type;
   template<typename A, typename B> view_type slice(A a, B b) const { return nothing();}
   nothing(...){}
   template<typename T1, typename T2, typename T3, typename T4> nothing(T1,T2,T3,T4){}
   template<typename T> void operator=(T) {}
   template<typename F> void assign_from_infty_evaluation(F f){}    
  };

  template<typename T, typename F> struct assign_to_F_of_infty;

 /*------------------------------------------------------------------------------
  * The implementation class for both the view and the regular gf class
  *-----------------------------------------------------------------------------*/
   template<typename MeshType, bool IsView> class gf_impl {

   public:

    typedef typename MeshType::gf_result_type value_element_type;
    typedef MeshType mesh_type;
    typedef MeshType domain_type;
    typedef std::vector<std::vector<std::string> > indices_type;

    typedef void has_view_type_tag; // so that triqs::lazy can put view in the expression trees. 
    typedef gf<MeshType,true>   view_type;
    typedef gf<MeshType,false>  non_view_type;
    typedef gf<MeshType,IsView> Derived;

    typedef typename mpl::if_c<MeshType::has_tail, gf<meshes::tail,false>, nothing>::type  tail_non_view_type; 
    typedef typename tail_non_view_type::view_type                                         tail_view_type; 
    typedef typename mpl::if_c<IsView, tail_view_type, tail_non_view_type>::type           tail_type;

    typedef tqa::array      <value_element_type,3,arrays::Option::Fortran>   data_non_view_type;
    typedef tqa::array_view <value_element_type,3,arrays::Option::Fortran>   data_view_type;
    typedef typename mpl::if_c<IsView, data_view_type, data_non_view_type>::type  data_type;

    typedef typename MeshType::index_type arg0_type;
    typedef tqa::matrix_view<value_element_type, arrays::Option::Fortran>       result_type;
    typedef tqa::matrix_view<const value_element_type, arrays::Option::Fortran> const_result_type;

   protected:
    mesh_type mesh;
    data_type data;
    tail_type tail;
   public:

    const indices_type indices;

    mesh_type const & the_mesh() const { return mesh;}

    data_view_type data_view() { return data;}
    const data_view_type data_view() const { return data;}

    tail_view_type tail_view() { return tail;}
    const tail_view_type tail_view() const { return tail;}

   protected:
    gf_impl (mesh_type const & mesh_, data_view_type const & data_, tail_view_type const & tail_, indices_type const & indices_) : 
     data(data_), mesh(mesh_), tail(tail_), indices(indices_){}

    //gf_impl(gf_impl const & x): mesh(x.mesh), data(x.data), tail(x.tail), indices(x.indices){}

    template<typename GfType> gf_impl(GfType const & x): mesh(x.the_mesh()), data(x.data_view()), tail(x.tail_view()), indices(x.indices){}
 
    gf_impl() {} 
 
    gf_impl (size_t N1, size_t N2, mesh_type const & mesh_, indices_type const & indices_) : 
      mesh(mesh_), data(N1,N2,mesh.len()), tail(N1,N2,mesh.mesh_tail,indices_), indices(indices_){}

   public:

    triqs::arrays::mini_vector<size_t,2> result_shape() const { return data.shape();}

    result_type       operator() ( arg0_type const & x)       { return data(range(),range(),x);}
    const_result_type operator() ( arg0_type const & x) const { return data(range(),range(),x);}

    tail_view_type       operator() ( domains::infty const & x)       { return tail;}
    const tail_view_type operator() ( domains::infty const & x) const { return tail;}

    TRIQS_LAZY_ADD_LAZY_CALL_WITH_VIEW(1,view_type);

    template<typename A, typename B> view_type slice(A a, B b) { 
     return view_type( gf_impl<MeshType,true>(mesh, data (range(a), range(b), range()), tail.slice(a,b), slice_indices(indices,range(a),range(b)) )); 
    }
    template<typename A, typename B> const view_type slice(A a, B b) const { 
     return view_type( gf_impl<MeshType,true>(mesh, data (range(a), range(b), range()), tail.slice(a,b), slice_indices(indices,range(a),range(b)) )); 
    }

   template<typename RHS> Derived & operator = (RHS const & rhs) { impl::assignment<Derived,RHS>::invoke(static_cast<Derived >(*this),rhs); return *this; } 

   // lazy_assignable // TO DO : the computation of the tail is not ready
   template<typename F> void set_from_function(F f) { 
    const size_t Nmax = this->data.shape()[2]; for (size_t u=0; u<Nmax; ++u) this->data(range(),range(),u) = f(u);
    compute_tail(tail,f); // the tail is only computed if it exists, or f maybe called for infty in undefined cases.
   }
   private:
   template<typename F> void compute_tail(nothing & t, F f) {}
   template<typename F> void compute_tail(gf_impl<meshes::tail, IsView> & t, F f) { std::cerr<<"impl tail assign"<<std::endl;} // t= f(domains::infty);
   public:

   /// Save the Green function in i omega_n (as 2 columns).
   void save(std::string file,  bool accumulate=false) const {}

   /// Load the GF
   void load(std::string file){}

   /// HDF5 saving ....

   private : 
   indices_type slice_indices(indices_type const & V, range const & R1, range const & R2){
    indices_type res(2); res[0].reserve(R1.last() - R1.first() + 1); res[1].reserve(R2.last() - R2.first() + 1); 
    for (size_t u =R1.first(); u<=R1.last(); ++u) res[0].push_back(V[0][u]);
    for (size_t u =R2.first(); u<=R2.last(); ++u) res[1].push_back(V[1][u]);
    return res;
   }
   friend class gf_impl<MeshType,!IsView>;
   };

 }// namespace impl 

 /**------------------------------------------------------------------------------
  * The View class of GF
  *-----------------------------------------------------------------------------*/
 template<typename MeshType> class gf<MeshType,true> : public impl::gf_impl<MeshType,true> {
  public :
   typedef impl::gf_impl<MeshType,true> base_type;

   // gf (mesh_type const & mesh_, data_view_type const & data_, tail_view_type const & tail_, indices_type const & indices_) : 
   // B(mesh_,data_,tail_,indices_) {}

   template<typename GfType> gf(GfType const & x): base_type(x) {};

   template<typename F> void set_from_function(F f) { base_type::set_from_function(f);} // bug of autodetection in triqs::lazy on gcc   
   
   std::ostream & print_for_lazy(std::ostream & out) const { return out<<"gf_view";}
 };

 /**------------------------------------------------------------------------------
  * The regular class of GF
  *-----------------------------------------------------------------------------*/
 template<typename MeshType> class gf<MeshType,false> : public impl::gf_impl<MeshType,false> {
  public : 
   typedef impl::gf_impl<MeshType,false> base_type;

   gf():base_type() {} 

   gf (size_t N1, size_t N2, typename base_type::mesh_type const & mesh_, typename base_type::indices_type const & indices_) : 
    base_type(N1,N2,mesh_,indices_) {}

   template<typename GfType> gf(GfType const & x): base_type(x) {} 
 };

}}}
#endif

