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
#include "tail.hpp"

namespace triqs { namespace gf { namespace local {
 
 template<typename MeshType> class gf;         // the value class
 template<typename MeshType> class gf_view;    // the view class
 template<typename E>        class gf_expr; // the proto expression of gf

 template<typename G> struct LocalGf : mpl::false_{};  // a boolean trait to identify the objects modelling the concept LocalGf
 template<typename M> struct LocalGf<gf<M> >     : mpl::true_{};
 template<typename M> struct LocalGf<gf_view<M> >: mpl::true_{};
 template<typename M> struct LocalGf<gf_expr<M> >: mpl::true_{};

 // ---------------------- implementation --------------------------------

 // traits to compute ...
 // special case : gf in time
 template<typename M> struct value_type_from_mesh_type { typedef dcomplex type;};
 template<typename M> struct tail_type_from_mesh_type { typedef tail type;};

 ///The regular class of GF
 template<typename MeshType,bool IsView> class gf_impl {
   public : 
   typedef MeshType mesh_type;
   typedef typename mesh_type::domain_type domain_type;
   static const bool has_tail = mesh_type::domain_type::has_tail;

   typedef typename value_type_from_mesh_type<mesh_type>::type value_type;
   typedef arrays::Option::Fortran storage_order;
   typedef arrays::array      <value_type,3,storage_order>                       data_non_view_type;
   typedef arrays::array_view <value_type,3,storage_order>                       data_view_type;
   typedef typename mpl::if_c<IsView, data_view_type, data_non_view_type>::type  data_type;

   typedef typename tail_type_from_mesh_type<mesh_type>::type                     tail_non_view_type; 
   typedef typename tail_non_view_type::view_type                                 tail_view_type; 
   typedef typename mpl::if_c<IsView, tail_view_type, tail_non_view_type>::type   tail_type;
  
   mesh_type const & mesh() const { return mymesh;}
   domain_type const & domain() const { return mymesh.domain();}
   //data_view_type data_view()             { return data;} 
   //const data_view_type data_view() const { return data;}

   typedef tqa::mini_vector<size_t,2> shape_type;
   shape_type shape() const { return shape_type(data.shape()[0], data.shape()[1]);} 

   typedef void has_view_type_tag; // so that triqs::lazy can put view in the expression trees. 
   typedef gf_view<MeshType> view_type;
   typedef gf<MeshType>      non_view_type;

  protected:
   mesh_type mymesh;
   data_type data;
   tail_type tail;
 
   gf_impl() {} // all arrays of zero size (empty)
   gf_impl(size_t N1, size_t N2, MeshType const & m, int tail_order_max): mymesh(m), data(N1,N2,m.size()), tail(N1,N2,-1,tail_order_max){ data()=0;}
   gf_impl (mesh_type const & m, data_view_type const & dat,tail_view_type const & t) : mymesh(m), data(dat), tail(t){}
   
   gf_impl(gf_impl const & x)                  : mymesh(x.mymesh), data(x.data), tail(x.tail){}// for clarity, would be synthetized anyway
   gf_impl(gf_impl<MeshType,!IsView> const & x): mymesh(x.mymesh), data(x.data), tail(x.tail){} 
   
   // orbital slice constructor. Only for view.
   template<bool V> gf_impl (gf_impl<MeshType,V> const & m, range r1, range r2) : 
     mymesh(m.mesh()),data(m.data(r1,r2, range())),tail(m.tail.slice(r1,r2)){ static_assert(IsView, "Internal Error");}
   
   // mesh slice constructor. Only for view.
   template<bool V> gf_impl (gf_impl<MeshType,V> const & m, typename MeshType::slice_arg_type arg) : 
     mymesh(m.mesh().slice(arg)),data(m.data),tail(m.tail){static_assert(IsView, "Internal Error");}

   friend class gf_impl<MeshType,!IsView>;

   // access to the data . Beware, we view it as a *matrix* NOT an array...
   tqa::matrix_view <value_type, storage_order> operator[](size_t u) { return data(range(),range(),u);}

   public:
   
   void operator = (gf_impl const & rhs) { mymesh = rhs.mymesh; data = rhs.data; tail = rhs.tail; }// for clarity, would be synthetized anyway
   void operator = (gf_impl<MeshType,!IsView> const & rhs) { mymesh = rhs.mymesh; data = rhs.data; tail = rhs.tail; }
   template<typename RHS> void operator = (RHS const & rhs) { // the general version  
    mymesh = rhs.mesh(); 
    shape_type sh = rhs.shape();
    tqa::resize_or_check_if_view( data, tqa::make_shape(sh[0],sh[1],mymesh.size()) );
    for (size_t u=0; u<mesh().size(); ++u)  (*this)[u] = rhs(mesh()[u]); 
    fill_tail1(mpl::bool_<has_tail>(), rhs); //tail = rhs (domains::infty());  
   }

   TRIQS_NVL_HAS_AUTO_ASSIGN(); 
   template<typename F> friend void triqs_nvl_auto_assign (gf_impl & x, F f) { // mesh is invariant in this case...
    for (size_t u=0; u<x.mesh().size(); ++u)  { x[u] = f(x.mesh()[u]); }
    fill_tail(mpl::bool_<gf_impl::has_tail>(), x.tail, f); 
    // if f is an expression, replace the placeholder with a simple tail. If f is a function callable on infty, 
    // it uses the fact that tail_non_view_type can be caster into domains::infty (See below).
   }

   private : // impl detail. We do not want that f or rhs is called when there is no tail...
   // using this overload, we guarantee that the call of f/rhs will not be *compiled* when not needed... 
   template<typename RHS> static void fill_tail ( mpl::false_, tail_type & t, RHS const & rhs) {} 
   template<typename RHS> static void fill_tail ( mpl::true_,  tail_type & t, RHS const & rhs) {t = rhs( tail::omega(t.shape(),t.size()));} 
   template<typename RHS> void fill_tail1 ( mpl::false_, RHS const & rhs) {} 
   template<typename RHS> void fill_tail1 ( mpl::true_,  RHS const & rhs) { this->tail = rhs (domains::infty()); }

   public : // call operators
   typedef arrays::matrix_view<value_type,       storage_order>  mv_type;
   typedef arrays::matrix_view<const value_type, storage_order>  const_mv_type;

   mv_type       operator() (meshes::mesh_pt<mesh_type> const & x)       { return data(range(),range(),x.i);}
   const_mv_type operator() (meshes::mesh_pt<mesh_type> const & x) const { return data(range(),range(),x.i);}

   tail_view_type       operator() ( domains::infty const & x)       {return tail;}
   const tail_view_type operator() ( domains::infty const & x) const {return tail;}

   typedef typename domain_type::point_type arg0_type;
   mv_type       operator() (arg0_type const & i)       { return this->mesh().interpolate(*this,i);}
   const_mv_type operator() (arg0_type const & i) const { return this->mesh().interpolate(*this,i);}

   TRIQS_LAZY_ADD_LAZY_CALL_WITH_VIEW(1,view_type);

   /// Slice in orbital space
   const view_type slice(range R1, range R2) const { return view_type(*this,R1,R2); } 

   /// Slice of the mesh : to be implemented
   view_type       slice_mesh(typename MeshType::slice_arg_type arg)       { return view_type (*this, arg);}
   view_type const slice_mesh(typename MeshType::slice_arg_type arg) const { return view_type (*this, arg);}

   /// Save the Green function in i omega_n (as 2 columns).
   void save(std::string file,  bool accumulate=false) const {}

   /// Load the GF
   void load(std::string file){}

   /// 
   friend void h5_write (tqa::h5::group_or_file fg, std::string subgroup_name, gf_impl const & g) {
    tqa::h5::group_or_file gr =  fg.create_group(subgroup_name);
    h5_write(gr,"data",g.data);
    h5_write(gr,"tail",g.tail);
    //h5_write(gr,"mesh",g.mymesh);
   }

   friend void h5_read  (tqa::h5::group_or_file fg, std::string subgroup_name, gf_impl & g){
    tqa::h5::group_or_file gr = fg.open_group(subgroup_name);
    h5_read(gr,"data",g.data);
    h5_read(gr,"tail",g.tail);
    //h5_read(gr,"mesh",g.mymesh);
   }

   //  BOOST Serialization
   friend class boost::serialization::access;
   template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
     ar & boost::serialization::make_nvp("data",data);
     ar & boost::serialization::make_nvp("tail",tail);
     ar & boost::serialization::make_nvp("mesh",mymesh);
    }

   /// print
   friend std::ostream & operator << (std::ostream & out, gf_impl const & x) { return out<<(IsView ? "gf_view": "gf");}
   friend std::ostream & triqs_nvl_formal_print(std::ostream & out, gf_impl const & x) { return out<<(IsView ? "gf_view": "gf");}
 };

 // ---------------------------------------------------------------------------------
 ///The regular class of GF
 template<typename MeshType> class gf : public gf_impl<MeshType,false> {
  typedef gf_impl<MeshType,false> B;
  public : 
  gf():B() {} 

  gf(gf const & g): B(g){}
  gf(gf_view<MeshType> const & g): B(g){} 
  template<typename GfType> gf(GfType const & x): B() { *this = x;} 

  gf(size_t N1, size_t N2, MeshType const & m, int tail_order_max = 3) : B(N1,N2,m,tail_order_max) {}

  using B::operator=;// or the default is = synthetized...
 };

 // ---------------------------------------------------------------------------------
 ///The View class of GF
 template<typename MeshType> class gf_view : public gf_impl<MeshType,true> {
  typedef gf_impl<MeshType,true> B;
  public :
  template<bool V> gf_view(gf_impl<MeshType,V> const & g): B(g){}
  template<bool V> gf_view(gf_impl<MeshType,V> const & g, range R1, range R2) : B(g,R1,R2){} // slice constructor 

  using B::operator=;// or the default is = synthetized...
 };

 // -------------------------------   Expression template for gf  --------------------------------------------------

 template< typename T> struct scalar_identification_trait : is_scalar_or_element<T> {}; // is_scalar_or_element defined in tail
 template< typename T> struct gf_identification_trait : LocalGf<T> {};
 static const int gf_arity=1;
 // this include in intentionnally IN the namespace
#include "../gf_proto.hpp"
}}}
#endif
