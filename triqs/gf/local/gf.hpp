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
 template<typename E>        class gf_expr;    // the proto expression of gf

 template<typename G> struct LocalGf : mpl::false_{};  // a boolean trait to identify the objects modelling the concept LocalGf
 template<typename M> struct LocalGf<gf<M> >     : mpl::true_{};
 template<typename M> struct LocalGf<gf_view<M> >: mpl::true_{};
 template<typename M> struct LocalGf<gf_expr<M> >: mpl::true_{};

 // ---------------------- implementation --------------------------------

 template<typename M> struct gf_value_type { typedef dcomplex type;};

 /// A common implementation class. Idiom : ValueView
 template<typename MeshType,bool IsView> class gf_impl {
   public : 
   typedef void has_view_type_tag;     // Idiom : ValueView  
   typedef gf_view<MeshType> view_type;
   typedef gf<MeshType>      non_view_type;
   
   typedef MeshType mesh_type;
   typedef typename mesh_type::domain_type domain_type;

   typedef typename gf_value_type<mesh_type>::type value_type;
   typedef arrays::Option::Fortran storage_order;
   typedef arrays::array      <value_type,3,storage_order>                       data_non_view_type;
   typedef arrays::array_view <value_type,3,storage_order>                       data_view_type;
   typedef typename mpl::if_c<IsView, data_view_type, data_non_view_type>::type  data_type;
   
   typedef typename mpl::if_c<IsView, tail_view, tail>::type   tail_type;
  
   mesh_type const & mesh() const { return _mesh;}
   domain_type const & domain() const { return _mesh.domain();}
   data_view_type data_view()             { return data;} 
   const data_view_type data_view() const { return data;}

   typedef tqa::mini_vector<size_t,2> shape_type;
   shape_type shape() const { return shape_type(data.shape()[0], data.shape()[1]);} 

  protected:
   mesh_type _mesh;
   data_type data;
   tail_type _tail;
 
   gf_impl() {} // all arrays of zero size (empty)
   gf_impl(size_t N1, size_t N2, MeshType const & m, int tail_order_max): _mesh(m), data(N1,N2,m.size()), _tail(N1,N2,-1,tail_order_max +2){ data()=0;}
   gf_impl(shape_type sh, MeshType const & m, tail_view const & t): _mesh(m), data(sh[0],sh[1],m.size()), _tail(t){ data()=0;}
   gf_impl (mesh_type const & m, data_view_type const & dat,tail_view const & t) : _mesh(m), data(dat), _tail(t){}
   
   gf_impl(gf_impl const & x)                  : _mesh(x._mesh), data(x.data), _tail(x._tail){}
   gf_impl(gf_impl<MeshType,!IsView> const & x): _mesh(x._mesh), data(x.data), _tail(x._tail){} 

   friend class gf_impl<MeshType,!IsView>;

   // access to the data . Beware, we view it as a *matrix* NOT an array... (crucial for assignment to scalars !)
   tqa::matrix_view <value_type, storage_order> _at (size_t u) { return data(tqa::range(),tqa::range(),u);}

   public:
   
   void operator = (gf_impl const & rhs) { _mesh = rhs._mesh; data = rhs.data; _tail = rhs._tail;}
   void operator = (gf_impl<MeshType,!IsView> const & rhs) { _mesh = rhs._mesh; data = rhs.data; _tail = rhs._tail;}
   template<typename RHS> void operator = (RHS const & rhs) { // the general version  
    _mesh = rhs.mesh(); 
    shape_type sh = rhs.shape();
    tqa::resize_or_check_if_view( data, tqa::make_shape(sh[0],sh[1],_mesh.size()) );
    for (size_t u=0; u<mesh().size(); ++u)  _at(u) = rhs(mesh()[u]); 
    _tail = rhs (domains::freq_infty());  
   }

   TRIQS_NVL_HAS_AUTO_ASSIGN(); 
   template<typename F> friend void triqs_nvl_auto_assign (gf_impl & x, F f) { // mesh is invariant in this case...
    for (size_t u=0; u<x.mesh().size(); ++u)  { x._at(u) = f(x.mesh()[u]); }
    x._tail = f( tail::omega(x._tail.shape(),x._tail.size()));
    // if f is an expression, replace the placeholder with a simple tail. If f is a function callable on freq_infty, 
    // it uses the fact that tail_non_view_type can be casted into domains::freq_infty 
   }
   
   // call operators
   typedef arrays::matrix_view<value_type,       storage_order>  mv_type;
   typedef arrays::matrix_view<const value_type, storage_order>  const_mv_type;

   mv_type       operator() (meshes::mesh_pt<mesh_type> const & x)       { return data(tqa::range(),tqa::range(),x.index);}
   const_mv_type operator() (meshes::mesh_pt<mesh_type> const & x) const { return data(tqa::range(),tqa::range(),x.index);}

   tail_view       operator() ( domains::freq_infty const & x)       {return _tail;}
   const tail_view operator() ( domains::freq_infty const & x) const {return _tail;}

   typedef typename domain_type::point_type arg0_type;
   mv_type       operator() (arg0_type const & i)       { return this->mesh().interpolate(*this,i);}
   const_mv_type operator() (arg0_type const & i) const { return this->mesh().interpolate(*this,i);}

   TRIQS_LAZY_ADD_LAZY_CALL_WITH_VIEW(1,view_type);

   /// Save the Green function in i omega_n (as 2 columns).
   void save(std::string file,  bool accumulate=false) const {}

   /// Load the GF
   void load(std::string file){}

   /// 
   friend void h5_write (tqa::h5::group_or_file fg, std::string subgroup_name, gf_impl const & g) {
    tqa::h5::group_or_file gr =  fg.create_group(subgroup_name);
    // Add the attribute
    h5_write(gr,"data",g.data);
    h5_write(gr,"tail",g._tail);
    //h5_write(gr,"mesh",g._mesh);
   }

   friend void h5_read  (tqa::h5::group_or_file fg, std::string subgroup_name, gf_impl & g){
    tqa::h5::group_or_file gr = fg.open_group(subgroup_name);
    // Check the attribute or throw 
    h5_read(gr,"data",g.data);
    h5_read(gr,"tail",g._tail);
    //h5_read(gr,"mesh",g._mesh);
   }

   //  BOOST Serialization
   friend class boost::serialization::access;
   template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
     ar & boost::serialization::make_nvp("data",data);
     ar & boost::serialization::make_nvp("tail",_tail);
     ar & boost::serialization::make_nvp("mesh",_mesh);
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
  gf(typename B::shape_type sh, MeshType const & m, tail_view const & t) : B(sh,m,t) {}
  using B::operator=;// or the default is = synthetized...
 };

 // ---------------------------------------------------------------------------------
 ///The View class of GF
 template<typename MeshType> class gf_view : public gf_impl<MeshType,true> {
  typedef gf_impl<MeshType,true> B;
  public :
  template<bool V> gf_view(gf_impl<MeshType,V> const & g): B(g){}
  template<typename D, typename T> gf_view (MeshType const & m, D const & dat,T const & t) : B(m,dat,t) {}
  using B::operator=;// or the default is = synthetized...
 };
   
 template<typename MeshType, bool V> 
  gf_view<MeshType> slice (gf_impl<MeshType,V> const & g, tqa::range r1, tqa::range r2) {
  return gf_view<MeshType>(g.mesh(), g.data_view()(r1,r2, tqa::range()), slice(g(domains::freq_infty()),r1,r2)); 
 }

 /*
 template<typename MeshType, bool V> slice (gf_impl<MeshType,V> const & g, typename MeshType::slice_arg_type arg) {
  return gf_view<MeshType>(slice_mesh(g.mesh()), g.data(r1,r2, tqa::range(???)), m._tail); 
 }
*/

 // -------------------------------   Expression template for gf  --------------------------------------------------

 template< typename T> struct scalar_identification_trait : is_scalar_or_element<T> {}; // is_scalar_or_element defined in tail
 template< typename T> struct gf_identification_trait : LocalGf<T> {};
 static const int gf_arity=1;
 // this include in intentionnally IN the namespace
#include "../gf_proto.hpp"
}}}
#endif
