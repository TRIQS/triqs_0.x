/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012-2013 by M. Ferrero, O. Parcollet
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
#ifndef TRIQS_GF_GFBASECLASS_H
#define TRIQS_GF_GFBASECLASS_H
#include <triqs/utility/first_include.hpp>
#include <triqs/utility/std_vector_expr_template.hpp>
#include <triqs/utility/factory.hpp>
#include "./tools.hpp"
#include <triqs/arrays/h5.hpp>
#include <vector>

namespace triqs { namespace gf { 
 using utility::factory;

 template<typename Descriptor> class gf;         // the value class
 template<typename Descriptor> class gf_view;    // the view class

 // The trait that "marks" the Green function
 TRIQS_DEFINE_CONCEPT_AND_ASSOCIATED_TRAIT(ImmutableGreenFunction);

 // ---------------------- implementation --------------------------------

 // we start with a few struct that helps to specialize some part of code
 // when the storage is arrays or is vector, depending on the descriptor.

 // indirection for gf assignment. Used later. Need it because the storage can be an array or a vector.
 template <typename G, typename Enable = void> struct gf_assigner { // vector case 
  template<typename RHS> 
   static void invoke(G & g, RHS && rhs) { 
   auto r = make_vector(rhs.data_view()); for (size_t i =0; i<g.data_view().size(); ++i) g.data_view()[i] = r[i];
   //g.data_view() = utility::factory<typename G::data_t>( utility::factory<typename G::data_non_view_t>(make_vector(rhs.data_view())));
  }
 };
 template <typename G> struct gf_assigner<G,ENABLE_IFC( tqa::is_array_or_view<typename G::data_t>::value)> { // array case
  template<typename RHS> static void invoke(G & g, RHS && rhs) { g.data_view() = rhs.data_view();}
 };
 template <typename G, typename Enable = void> struct gf_assigner_to_scalar { // vector case. 
  template<typename RHS> static void invoke(G & g, RHS const & rhs) { for (size_t i =0; i<g.data_view().size(); ++i) g.data_view()[i] = rhs; }
 };
 template <typename G> struct gf_assigner_to_scalar<G,ENABLE_IFC( tqa::is_array_or_view<typename G::data_t>::value)> { //array
  template<typename RHS> static void invoke(G & g, RHS const & rhs) { g.data_view() = rhs;}
 };

 //------------------------------------------------------------------------

 // indirection for calling the h5 read/write for the storage
 // for e.g. block, it is handled by the descriptor itself.
 // otherwise it is the default...
 template<typename D, typename Enable=void> struct gf_h5_caller {
  static void write(h5::group g, std::string const & s, typename D::mesh_t const & mesh, typename D::storage_t const & data) { 
   h5_write(g,"data",data);
  }
  static void read(h5::group g, std::string const & s, typename D::mesh_t const & mesh, typename D::storage_t & data) { 
   h5_read(g,"data",data);
  }
 };

 // detect if D has the h5_data_read function
 template<typename D>
  struct gf_h5_caller<D, typename D::h5_use_special_read_write>  {
   static void write(h5::group g, std::string const & s, typename D::mesh_t const & mesh, typename D::storage_t const & data) {
    D::h5_data_write(g,s,mesh,data);
   }

   static void read(h5::group g, std::string const & s, typename D::mesh_t const & mesh, typename D::storage_t & data) { 
    D::h5_data_read(g,s,mesh,data);
   }
  };

 //------------------------------------------------------------------------

 // a little piece of code to get the value of the gf, from a vector, an array 1d, or 3d
 // depending of the descriptor. The call to () is [] is different, so need to specialize 
 template<typename T, typename Enable = void, typename Enable2 = void> struct gf_data_getter;
 template<typename T> struct gf_data_getter<T, ENABLE_IFC( tqa::is_array_or_view<T>::value), ENABLE_IFC(T::rank!=1)> { 
  static auto invoke(T & X, size_t i)       -> decltype(X(tqa::ellipsis(),i)) { return X(tqa::ellipsis(),i);}
  static auto invoke(T const & X, size_t i) -> decltype(X(tqa::ellipsis(),i)) { return X(tqa::ellipsis(),i);}
 };
 template<typename T> struct gf_data_getter<T, ENABLE_IFC( tqa::is_array_or_view<T>::value), ENABLE_IFC(T::rank==1)> { 
  static auto invoke(T & X, size_t i)       -> decltype(X(i)) { return X(i);}
  static auto invoke(T const & X, size_t i) -> decltype(X(i)) { return X(i);}
 };
 template<typename T> struct gf_data_getter<std::vector<T>, void, void> { 
  static auto invoke(std::vector<T> &X, size_t i)       -> decltype(X[i]) { return X[i];}
  static auto invoke(std::vector<T> const &X, size_t i) -> decltype(X[i]) { return X[i];}
 };

 //------------------------------------------------------------------------
// a little function to specialize the rebind a vector of views ...
 
 template<typename T> void gf_data_rebinder(T & a, T const & rhs) { a.rebind(rhs);}
 template<typename T> void gf_data_rebinder(std::vector<T> & a, std::vector<T> const & rhs){a.clear(); for (auto & x : rhs) a.push_back(x);}

 //------------------------------------------------------------------------
 /// A common implementation class. Idiom : ValueView
 template<typename Descriptor,bool IsView> class gf_impl : TRIQS_MODEL_CONCEPT(ImmutableGreenFunction), Descriptor::tag {
  public : 

   typedef Descriptor                             descriptor_t;
   typedef typename Descriptor::mesh_t            mesh_t;
   typedef typename mesh_t::domain_t              domain_t;
   typedef typename Descriptor::target_t          target_t;
   typedef typename target_t::value_type          value_t;
   typedef typename target_t::value_type          value_type; // CLEAN THIS EVERYWHERE !!!

   typedef typename mesh_t::mesh_point_t          mesh_point_t;
   typedef typename mesh_t::index_t               mesh_index_t;
   typedef typename Descriptor::singularity_t     singularity_non_view_t;
   typedef typename Descriptor::symmetry_t        symmetry_t;
   typedef typename Descriptor::indices_t         indices_t;

   typedef typename view_type_if_exists_else_type<target_t>::type          target_view_t;
   typedef typename view_type_if_exists_else_type<singularity_non_view_t>::type     singularity_view_t;

   typedef void has_view_type_tag;     // Idiom : ValueView  
   typedef gf_view<Descriptor> view_type;
   typedef gf<Descriptor>      non_view_type;

   typedef typename Descriptor::storage_t                                  data_non_view_t;
   typedef typename Descriptor::storage_view_t                             data_view_t;
   typedef typename mpl::if_c<IsView, data_view_t, data_non_view_t>::type  data_t;

   typedef typename mpl::if_c<IsView, singularity_view_t, singularity_non_view_t>::type   singularity_t;

   mesh_t const & mesh() const            { return _mesh;}
   domain_t const & domain() const        { return _mesh.domain();}
   data_t &  data_view()                  { return data;}
   data_t const & data_view() const       { return data;}
   singularity_view_t singularity_view()  { return singularity;} 
   const singularity_view_t singularity_view() const { return singularity;}

   symmetry_t const & symmetry() const { return _symmetry;}
   indices_t const & indices() const   { return _indices;}

  protected:
   mesh_t _mesh;
   data_t data;
   singularity_t singularity;
   symmetry_t _symmetry;
   indices_t _indices;
   typename Descriptor::evaluator _evaluator;
   typename Descriptor::bracket_evaluator _bracket_evaluator;

   // Constructors : protected, see gf/gf_view later for public one
   gf_impl() {} // all arrays of zero size (empty)
  public : //everyone can make a copy (for clef lib in particular, this is needed)
   gf_impl(gf_impl const & x) = default;
   // do the move constructor
  protected:
   gf_impl(gf_impl<Descriptor,!IsView> const & x): _mesh(x.mesh()), data(utility::factory<data_t>(x.data_view())), 
   singularity(x.singularity_view()), _symmetry(x.symmetry()), _indices(x.indices()){} 

   // from the data directly
   gf_impl(mesh_t const & m, data_view_t const & dat, singularity_view_t const & ad, symmetry_t const & s, indices_t const & ind ) : 
    _mesh(m), data(dat), singularity(ad),_symmetry(s), _indices(ind){}

   gf_impl(mesh_t const & m, data_t && dat, singularity_view_t const & ad, symmetry_t const & s, indices_t const & ind ) : 
    _mesh(m), data(std::move(dat)), singularity(ad),_symmetry(s), _indices(ind) {}

   void operator = (gf_impl const & rhs) = delete; // done in derived class.
   //bool is_empty() const { return data.is_empty();} // FIX THIS : does not work for vector !

  public:

   // ------------- All the call operators -----------------------------

   // First, a simple () returns a view, like for an array...
   view_type operator()() const { return *this;}

   /// Calls are (perfectly) forwarded to the Descriptor::operator(), except mesh_point_t and when 
   /// there is at least one lazy argument ...
   template<typename Arg0, typename... Args >    // match any argument list, picking out the first type : () is not permitted
    typename std::add_const<
    typename boost::lazy_disable_if<  // disable the template if one the following conditions it true 
    boost::mpl::or_< // starting condition [OR] 
    boost::is_base_of< typename std::remove_reference<Arg0>::type, mesh_point_t>,  // Arg0 is (a & or a &&) to a mesh_point_t 
    clef::is_any_lazy<Arg0, Args...>                          // One of Args is a lazy expression
     >,                                                       // end of OR 
    std::result_of<typename Descriptor::evaluator(mesh_t, data_t, singularity_t, Arg0, Args...)> // what is the result type of call
     >::type     // end of lazy_disable_if
     >::type // end of add_Const 
     operator() (Arg0&& arg0, Args&&... args) const { return _evaluator(_mesh, data, singularity, std::forward<Arg0>( arg0), std::forward<Args>(args)...); }

   // Interaction with the CLEF library : calling the gf with any clef expression as argument build a new clef expression
   template<typename Arg0, typename ...Args>
    //typename boost::lazy_enable_if<    // enable the template if
    //clef::is_any_lazy<Arg0, Args...>,  // One of Args is a lazy expression
    typename clef::result_of::make_expr_call<view_type,Arg0, Args...>::type
    //typename clef::result_of::make_expr_call<gf_impl,Arg0, Args...>::type
    //clef::result_of::make_expr_call<view_type,Args...>
    //>::type     // end of lazy_enable_if 
    operator()(Arg0 arg0, Args... args) const { 
     static_assert(sizeof...(Args) == Descriptor::arity-1, "Incorrect number of variable in call");
     return clef::make_expr_call(view_type(*this),arg0, args...);
    }

   /// A direct access to the grid point 
   target_view_t operator() (mesh_point_t const & x)       { return gf_data_getter<data_t>::invoke(data, x.m->index_to_linear(x.index));}

   /// A direct access to the grid point (const version)
   target_view_t operator() (mesh_point_t const & x) const { return gf_data_getter<data_t>::invoke(data, x.m->index_to_linear(x.index));}

   /// A direct access to the grid point 
   template<typename... Args> 
    target_view_t on_mesh (Args&&... args)       { return gf_data_getter<data_t>::invoke(data,_mesh.index_to_linear(mesh_index_t(std::forward<Args>(args)...)));}

   /// A direct access to the grid point (const version)
   template<typename... Args> 
    const target_view_t on_mesh (Args&&... args) const { return gf_data_getter<data_t>::invoke(data, _mesh.index_to_linear(mesh_index_t(std::forward<Args>(args)...)));}

  private:
   struct _on_mesh_wrapper {
    gf_impl const & f; _on_mesh_wrapper (gf_impl const & _f) : f(_f) {}
    template <typename... Args> target_view_t operator ()(Args && ... args) const { return f.on_mesh(std::forward<Args>(args)...);} 
   };
   _on_mesh_wrapper friend on_mesh(gf_impl const & f) { return f;}

  public:

   /// [] Calls are (perfectly) forwarded to the Descriptor::operator[]
   //except when there is at least one lazy argument ...
   template<typename Arg >   
    typename boost::lazy_disable_if<  // disable the template if one the following conditions it true 
    clef::is_any_lazy<Arg>,                          // One of Args is a lazy expression
    std::result_of<typename Descriptor::bracket_evaluator(mesh_t, data_t &, singularity_t &, Arg)> // what is the result type of call
     >::type     // end of lazy_disable_if 
     operator[] (Arg&& arg) {return _bracket_evaluator(_mesh, this->data, singularity, std::forward<Arg>( arg));}

   /// [] Calls are (perfectly) forwarded to the Descriptor::operator[]
   //except when there is at least one lazy argument ...
   template<typename Arg >   
    typename boost::lazy_disable_if<  // disable the template if one the following conditions it true 
    clef::is_any_lazy<Arg>,                          // One of Args is a lazy expression
    std::result_of<typename Descriptor::bracket_evaluator(mesh_t, data_t const &, singularity_t&, Arg)> // what is the result type of call
     >::type     // end of lazy_disable_if 
     operator[] (Arg&& arg) const {return _bracket_evaluator(_mesh, data, singularity, std::forward<Arg>( arg));}

   // Interaction with the CLEF library : calling the gf with any clef expression as argument build a new clef expression
   template<typename Arg>
    typename boost::lazy_enable_if<    // enable the template if
    clef::is_any_lazy<Arg>,  // One of Args is a lazy expression
    clef::result_of::make_expr_call<view_type,Arg>
     >::type     // end of lazy_enable_if 
     operator[](Arg const & arg) const { return  clef::make_expr_call(view_type(*this),arg);} //clef::lazy_call<view_type>(*this)(arg); 

   /// Write into HDF5
   friend void h5_write (h5::group fg, std::string subgroup_name, gf_impl const & g) {
    auto gr =  fg.create_group(subgroup_name);
    // Add the attribute
    gf_h5_caller<Descriptor>::write(gr,"data",g.mesh(), g.data);
    h5_write(gr,"singularity",g.singularity);
    h5_write(gr,"mesh",g._mesh);
    h5_write(gr,"symmetry",g._symmetry);
    h5_write(gr,"indices",g._indices);
   }

   /// Read from HDF5
   friend void h5_read  (h5::group fg, std::string subgroup_name, gf_impl & g){
    auto gr = fg.open_group(subgroup_name);
    // Check the attribute or throw 
    gf_h5_caller<Descriptor>::read(gr,"data",g.mesh(),g.data);
    h5_read(gr,"singularity",g.singularity);
    h5_read(gr,"mesh",g._mesh);
    h5_read(gr,"symmetry",g._symmetry);
    h5_read(gr,"indices",g._indices);
   }

   //  BOOST Serialization
   friend class boost::serialization::access;
   template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
     ar & boost::serialization::make_nvp("data",data);
     ar & boost::serialization::make_nvp("singularity",singularity);
     ar & boost::serialization::make_nvp("mesh",_mesh);
     ar & boost::serialization::make_nvp("symmetry",_symmetry);
     ar & boost::serialization::make_nvp("indices",_indices);
    }

   /// print
   friend std::ostream & operator << (std::ostream & out, gf_impl const & x) { return out<<(IsView ? "gf_view": "gf");}
   friend std::ostream & triqs_nvl_formal_print(std::ostream & out, gf_impl const & x) { return out<<(IsView ? "gf_view": "gf");}
 };

 // ---------------------------------------------------------------------------------
 ///The regular class of GF
 template<typename Descriptor> class gf : public gf_impl<Descriptor,false> {
  typedef gf_impl<Descriptor,false> B;
  public : 

  gf():B() {} 
  gf(gf const & g): B(g){}
  gf(gf && g): B(std::move(g)){}
  gf(gf_view<Descriptor> const & g): B(g){} 
  template<typename GfType> gf(GfType const & x): B() { *this = x;} 
  
  template<typename DATA_VIEW_TYPE> // anything from which the factory can make the data ...
   gf(typename B::mesh_t const & m, 
     DATA_VIEW_TYPE const & dat, 
     //typename B::data_view_t const & dat, 
     typename B::singularity_view_t const & si, 
     typename B::symmetry_t const & s , 
     typename B::indices_t const & ind = typename B::indices_t () ) : 
    B(m,factory<typename B::data_t>(dat),si,s, ind) {}

  void operator = (gf const & rhs) { *this = gf(rhs);} // use move =

  friend void swap ( gf & a, gf & b) { 
   swap(a._mesh, b._mesh); swap(a.data, b.data); swap (a.singularity,b.singularity); swap(a._symmetry,b._symmetry);swap(a._indices,b._indices);
  }

  void operator = (gf && rhs) { //swap(*this, rhs);
   this->_mesh = rhs._mesh; 
   this->data = rhs.data;
   this->singularity = rhs.singularity;
   this->_symmetry = rhs._symmetry; this->_indices = rhs._indices;
  }

  template<typename RHS> void operator = (RHS && rhs) { 
   this->_mesh = rhs.mesh(); 
   gf_assigner<gf>::invoke(*this, std::forward<RHS>(rhs));
   //this->data = rhs.data_view(); 
   this->singularity = rhs.singularity_view();
   //   this->_symmetry = rhs.symmetry(); this->_indices = rhs._indices();
  }
 };

 // ---------------------------------------------------------------------------------
 ///The View class of GF
 template<typename Descriptor> class gf_view : public gf_impl<Descriptor,true> {
  typedef gf_impl<Descriptor,true> B;

  public :

  void rebind( gf_view const &X) { 
   this->_mesh = X._mesh; this->_symmetry = X._symmetry; this->_indices = X._indices; 
   gf_data_rebinder(this->data,X.data); 
   //this->data.rebind(X.data); 
   this->singularity.rebind(X.singularity); 
  }

  gf_view(gf_view const & g): B(g){}
  template<bool V> gf_view(gf_impl<Descriptor,V> const & g): B(g){}

  template<typename D, typename T> 
   gf_view (typename B::mesh_t const & m, 
     D const & dat,T const & t,typename B::symmetry_t const & s, 
     typename B::indices_t const & ind = typename B::indices_t () ) : B(m,factory<typename B::data_t>(dat),t,s,ind) {}

  void operator = (gf_view const & rhs)  { triqs_gf_view_assign_delegation(*this,rhs);}

  template<typename RHS> void operator = (RHS const & rhs) { triqs_gf_view_assign_delegation(*this,rhs);}

  // Interaction with the CLEF library : auto assignment of the gf (gf(om_) << expression fills the functions by evaluation of expression)
  template<typename F> friend void triqs_clef_auto_assign (gf_view & x, F f) { Descriptor::assign_from_expression(x._mesh,x.data, x.singularity,f);}
 };

 // delegate = so that I can overload it for specific RHS...
 template<typename Descriptor, typename RHS> 
  void triqs_gf_view_assign_delegation( gf_view<Descriptor> & g, RHS && rhs) { 
   if (!(g.mesh()  == rhs.mesh()))  TRIQS_RUNTIME_ERROR<<"Gf Assignment in View : incompatible mesh";
   //if (!(g.shape() == rhs.shape())) TRIQS_RUNTIME_ERROR<<"Gf Assignment in View : incompatible target shape";
   gf_assigner<gf_view<Descriptor>>::invoke(g, rhs);
   //g.data_view() = rhs.data_view();
   g.singularity_view() = rhs.singularity_view();
  }

 template<typename Descriptor, typename T>
  ENABLE_IF(arrays::is_scalar<T>) triqs_gf_view_assign_delegation( gf_view<Descriptor> & g, T const & x) {
   gf_assigner_to_scalar<gf_view<Descriptor>>::invoke(g ,x);
   //g.data_view() = x; 
   g.singularity_view() = x;
  }

 // tool for lazy transformation
 template<typename Tag, typename D> struct gf_keeper{ gf_view<D> g; gf_keeper (gf_view<D> const & g_) : g(g_) {} };

 // ---------------------------------- slicing ------------------------------------

 template<typename Descriptor, bool V, typename... Args> 
  gf_view<Descriptor> slice_target (gf_impl<Descriptor,V> const & g, Args... args) {
   return gf_view<Descriptor>(g.mesh(), g.data_view()(args... , tqa::range()), slice_target (g.singularity_view(),args...), g.symmetry()); 
  }

 template<typename Descriptor, bool V, typename... Args> 
  gf_view<Descriptor> slice_mesh (gf_impl<Descriptor,V> const & g, Args... args) {
   return gf_view<Descriptor>(g.mesh().slice(args...), g.data_view()(arrays::ellipsis(), g.mesh().slice_get_range(args...)), g.singularity_view(), g.symmetry()); 
  }

}}

#include "./gf_expr.hpp"
#endif
