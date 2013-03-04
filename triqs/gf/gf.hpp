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
#include <triqs/arrays/view_proxy.hpp>
#include <vector>

namespace triqs { namespace gf {
 using utility::factory;

 template<typename Descriptor> class gf;         // the value class
 template<typename Descriptor> class gf_view;    // the view class

 template <typename Descriptor> struct gf_factories;
 template <typename Descriptor, typename ... U> gf<Descriptor> make_gf(U && ... x) { return gf_factories<Descriptor>::make_gf(std::forward<U>(x)...);}
 template <typename Descriptor, typename ... U> gf_view<Descriptor> make_gf_view(U && ... x) { return gf_factories<Descriptor>::make_gf_view(std::forward<U>(x)...);}

 template<typename Descriptor, typename G> struct evaluator; // Dispatch of all the () call except mesh_point
 template<typename Descriptor> struct data_proxy;            // 

 // The trait that "marks" the Green function
 TRIQS_DEFINE_CONCEPT_AND_ASSOCIATED_TRAIT(ImmutableGreenFunction);

 // ---------------------- implementation --------------------------------

 // detects if the descriptor has a special h5 read/write method
 template<typename D, typename Enable=void> struct has_special_h5_read_write                    : std::false_type {};
 template<typename D> struct has_special_h5_read_write<D, typename D::has_special_h5_read_write_tag>: std::true_type {};

 //------------------------------------------------------------------------

 // a little piece of code to get the value of the gf, from a vector, an array 1d, or 3d
 // depending of the descriptor. The call to () is [] is different, so need to specialize
 template<typename T, typename Enable = void> struct gf_data_getter;

 template<typename T> struct gf_data_getter<T, ENABLE_IFC( tqa::is_array_or_view<T>::value && (T::rank!=1))> {
  typedef arrays::view_proxy<T,2> r_type;
  typedef arrays::const_view_proxy<T,2> cr_type;
  r_type  operator()(T       & data, size_t i)       { return r_type(data,i); } 
  cr_type operator()(T const & data, size_t i) const { return cr_type(data,i); } 
 };

 template<typename T> struct gf_data_getter<T, ENABLE_IFC( tqa::is_array_or_view<T>::value && (T::rank==1))> {
  //typedef typename T::value_type & r_type;
  //typedef typename T::value_type const & cr_type;
  auto operator()(T       & data,size_t i)       -> decltype(data(i)) { return data(i);}
  auto operator()(T const & data,size_t i) const -> decltype(data(i)) { return data(i);}
 };

 template<typename T> struct gf_data_getter<std::vector<T>, void> {
  //typedef T       & r_type;
  //typedef T const & cr_type;
  T       &  operator()(std::vector<T> &       data, size_t i)  { return data[i];}
  T const &  operator()(std::vector<T> const & data, size_t i)  { return data[i];}
 };

 //------------------------------------------------------------------------
 /// A common implementation class. Idiom : ValueView
 template<typename Descriptor,bool IsView> class gf_impl : TRIQS_MODEL_CONCEPT(ImmutableGreenFunction), Descriptor::tag {
  public :

   typedef void has_view_type_tag;     // Idiom : ValueView
   typedef gf_view<Descriptor> view_type;
   typedef gf<Descriptor>      non_view_type;

   typedef Descriptor                             descriptor_t;
   typedef typename Descriptor::mesh_t            mesh_t;
   typedef typename mesh_t::domain_t              domain_t;
   typedef typename mesh_t::mesh_point_t          mesh_point_t;
   typedef typename mesh_t::index_t               mesh_index_t;
   typedef typename Descriptor::symmetry_t        symmetry_t;
   typedef typename Descriptor::indices_t         indices_t;

   typedef typename Descriptor::storage_t                                  data_non_view_t;
   typedef typename Descriptor::storage_view_t                             data_view_t;
   typedef typename std::conditional<IsView, data_view_t, data_non_view_t>::type  data_t;

   typedef typename Descriptor::singularity_t     singularity_non_view_t;
   typedef typename view_type_if_exists_else_type<singularity_non_view_t>::type     singularity_view_t;
   typedef typename std::conditional<IsView, singularity_view_t, singularity_non_view_t>::type   singularity_t;

   mesh_t const & mesh() const            { return _mesh;}
   domain_t const & domain() const        { return _mesh.domain();}
   data_t &  data_view()                  { return data;}
   data_t const & data_view() const       { return data;}
   singularity_view_t singularity_view()  { return singularity;}
   const singularity_view_t singularity_view() const { return singularity;}

   symmetry_t const & symmetry() const { return _symmetry;}
   indices_t const & indices() const { return _indices;}

  protected:
   mesh_t _mesh;
   data_t data;
   singularity_t singularity;
   symmetry_t _symmetry;
   indices_t _indices;
   typedef evaluator<Descriptor,gf_impl> evaluator_t;
   evaluator_t evaluator_;
   typedef gf_data_getter<data_t> data_getter_t;
   data_getter_t data_getter;
  public:
   std::integral_constant<bool, arrays::is_array_or_view<data_t>::value> storage_is_array;
  protected:

   // --------------------------------Constructors -----------------------------------------------
   // protected, see gf/gf_view later for public one
   gf_impl() : evaluator_(*this) {} // all arrays of zero size (empty) 
  public : //everyone can make a copy (for clef lib in particular, this is needed)
   gf_impl(gf_impl const & x) : _mesh(x.mesh()), data(factory<data_t>(x.data_view())),
   singularity(x.singularity_view()), _symmetry(x.symmetry()), _indices(x.indices()),evaluator_(*this){}

   gf_impl(gf_impl && ) = default;
  protected:
   gf_impl(gf_impl<Descriptor,!IsView> const & x): _mesh(x.mesh()), data(factory<data_t>(x.data_view())),
   singularity(x.singularity_view()), _symmetry(x.symmetry()), _indices(x.indices()),evaluator_(*this){}

   // from the data directly
   gf_impl(mesh_t const & m, data_view_t const & dat, singularity_view_t const & ad, symmetry_t const & s, indices_t const & ind ) :
    _mesh(m), data(dat), singularity(ad),_symmetry(s), _indices(ind),evaluator_(*this){}

   gf_impl(mesh_t const & m, data_t && dat, singularity_view_t const & ad, symmetry_t const & s, indices_t const & ind ) :
    _mesh(m), data(std::move(dat)), singularity(ad),_symmetry(s), _indices(ind), evaluator_(*this) {}

   // -------------------------------- assigner -----------------------------------------------

   void operator = (gf_impl const & rhs) = delete; // done in derived class.

  public:
   // array is true, vector is false.
   // warning : dimension should match, no resize ??
   template<typename RHS> void _data_assigner_no_resize (RHS && rhs, std::true_type) { data() = rhs.data_view();}
   template<typename RHS> void _data_assigner_no_resize (RHS && rhs, std::false_type) {
    auto r = make_vector(rhs.data_view()); 
    if (data.size() !=r.size()) TRIQS_RUNTIME_ERROR << "Size mismatch in gf assignment";
    for (size_t i =0; i<data.size(); ++i) data[i] = r[i];
   }

   // data can change size, we need to reset the data_getter
   // This is the only point where the arrays can be reshaped, hence views been invalidated !! 
   template<typename RHS> void _data_assigner_with_resize( RHS && rhs, std::true_type){data = rhs.data_view();}
   template<typename RHS> void _data_assigner_with_resize( RHS && rhs, std::false_type){data = factory<data_t>(rhs.data_view());}

   template<typename RHS> void _data_assigner_to_scalar( RHS && rhs, std::true_type){ data() = rhs; }
   template<typename RHS> void _data_assigner_to_scalar( RHS && rhs, std::false_type){for (size_t i =0; i<data.size(); ++i) data[i] = rhs;}

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
    std::result_of<evaluator_t(Arg0, Args...)> // what is the result type of call
     >::type     // end of lazy_disable_if
     >::type // end of add_Const
     operator() (Arg0&& arg0, Args&&... args) const { return evaluator_(std::forward<Arg0>( arg0), std::forward<Args>(args)...); }
   //operator() (Arg0&& arg0, Args&&... args) const { return _evaluator(_mesh, data, singularity, std::forward<Arg0>( arg0), std::forward<Args>(args)...); }

   // Interaction with the CLEF library : calling the gf with any clef expression as argument build a new clef expression
   template<typename Arg0, typename ...Args>
    //typename boost::lazy_enable_if<    // enable the template if
    //clef::is_any_lazy<Arg0, Args...>,  // One of Args is a lazy expression
    typename clef::result_of::make_expr_call<view_type,Arg0, Args...>::type
    //typename clef::result_of::make_expr_call<gf_impl,Arg0, Args...>::type
    //clef::result_of::make_expr_call<view_type,Args...>
    //>::type     // end of lazy_enable_if
    operator()(Arg0 arg0, Args... args) const {
     static_assert(sizeof...(Args) == evaluator_t::arity-1, "Incorrect number of variable in call");
     return clef::make_expr_call(view_type(*this),arg0, args...);
    }

   /// A direct access to the grid point (const version)
   //auto operator() (mesh_point_t const & x) const DECL_AND_RETURN(data_getter(data, x.m->index_to_linear(x.index)));

   typedef typename std::result_of<data_getter_t(data_t       &,size_t)>::type r_type;
   typedef typename std::result_of<data_getter_t(data_t const &,size_t)>::type cr_type;

   cr_type operator() (mesh_point_t const & x) const { return data_getter(data, x.m->index_to_linear(x.index));}

   /// A direct access to the grid point
   template<typename... Args>
    r_type on_mesh (Args&&... args) { return data_getter(data,_mesh.index_to_linear(mesh_index_t(std::forward<Args>(args)...)));}

   /// A direct access to the grid point (const version)
   template<typename... Args>
    cr_type on_mesh (Args&&... args) const { return data_getter(data,_mesh.index_to_linear(mesh_index_t(std::forward<Args>(args)...)));}

  private:
   struct _on_mesh_wrapper {
    gf_impl const & f; _on_mesh_wrapper (gf_impl const & _f) : f(_f) {}
    template <typename... Args> typename data_getter_t::cr_type operator ()(Args && ... args) const { return f.on_mesh(std::forward<Args>(args)...);}
    template <typename... Args> typename data_getter_t::r_type  operator ()(Args && ... args)       { return f.on_mesh(std::forward<Args>(args)...);}
   };
   _on_mesh_wrapper friend on_mesh(gf_impl const & f) { return f;}

  public:
   r_type  operator[] ( mesh_index_t const & arg)       { return data_getter(data,_mesh.index_to_linear(arg));}
   cr_type operator[] ( mesh_index_t const & arg) const { return data_getter(data,_mesh.index_to_linear(arg));}

   /// A direct access to the grid point
   r_type  operator[] (mesh_point_t const & x) { return data_getter(data, _mesh.index_to_linear(x.index));}
   cr_type operator[] (mesh_point_t const & x) const { return data_getter(data, _mesh.index_to_linear(x.index));}

   // Interaction with the CLEF library : calling the gf with any clef expression as argument build a new clef expression
   template<typename Arg>
    typename boost::lazy_enable_if<    // enable the template if
    clef::is_any_lazy<Arg>,  // One of Args is a lazy expression
    clef::result_of::make_expr_call<view_type,Arg>
     >::type     // end of lazy_enable_if
     operator[](Arg const & arg) const { return clef::make_expr_call(view_type(*this),arg);}

   //----------------------------- HDF5 -----------------------------

  private : // indirection if descrptor has a special read write
   void __h5_write (h5::group g, std::string const & s, std::false_type) const { h5_write(g,"data",data); }
   void __h5_write (h5::group g, std::string const & s, std::true_type)  const { Descriptor::h5_data_write(g,s,mesh,data);}
   void __h5_read  (h5::group g, std::string const & s, std::false_type) { h5_read(g,"data",data); }
   void __h5_read  (h5::group g, std::string const & s, std::true_type)  { Descriptor::h5_data_read(g,s,mesh,data);}
  public:

   /// Write into HDF5
   friend void h5_write (h5::group fg, std::string subgroup_name, gf_impl const & g) {
    auto gr =  fg.create_group(subgroup_name);
    // Add the attribute
    g.__h5_write (gr,"data", has_special_h5_read_write<Descriptor>());
    h5_write(gr,"singularity",g.singularity);
    h5_write(gr,"mesh",g._mesh);
    h5_write(gr,"symmetry",g._symmetry);
    h5_write(gr,"indices",g._indices);
   }

   /// Read from HDF5
   friend void h5_read  (h5::group fg, std::string subgroup_name, gf_impl & g){
    auto gr = fg.open_group(subgroup_name);
    // Check the attribute or throw
    g.__h5_read (gr,"data", has_special_h5_read_write<Descriptor>());
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
     typename B::singularity_view_t const & si,
     typename B::symmetry_t const & s ,
     typename B::indices_t const & ind = typename B::indices_t () ) :
    B(m,factory<typename B::data_t>(dat),si,s, ind) {}

  friend void swap ( gf & a, gf & b) {
   std::swap(a._mesh, b._mesh); swap(a.data, b.data); std::swap (a.singularity,b.singularity); std::swap(a._symmetry,b._symmetry);
   std::swap(a._indices,b._indices); 
  }

  void operator = (gf const & rhs) { *this = gf(rhs);} // use move =
  void operator = (gf && rhs) { swap(*this, rhs); }

  template<typename RHS> void operator = (RHS && rhs) {
   this->_mesh = rhs.mesh();
   this->_data_assigner_with_resize(std::forward<RHS>(rhs), this->storage_is_array);
   //this->data = rhs.data_view();
   this->singularity = rhs.singularity_view();
   // to be implemented : there is none in the gf_expr in particular....
   //this->_symmetry = rhs.symmetry(); this->_indices = rhs._indices();
  }
 };

 // ---------------------------------------------------------------------------------
 ///The View class of GF
 template<typename Descriptor> class gf_view : public gf_impl<Descriptor,true> {
  typedef gf_impl<Descriptor,true> B;
  // a little function to specialize the rebind a vector of views (false) or arrays (true)
  void _rebinder( gf_view const &X, std::true_type)  { this->data.rebind(X.data); }
  void _rebinder( gf_view const &X, std::false_type) { this->data.clear(); for (auto & x : X.data) this->data.push_back(x);}
  public :

  void rebind( gf_view const &X) {
   this->_mesh = X._mesh; this->_symmetry = X._symmetry; this->_indices = X._indices;
   _rebinder(X, this->storage_is_array); // need to change the trait to std::integrall...
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
   g._data_assigner_no_resize(std::forward<RHS>(rhs), g.storage_is_array); //g.data_view() = rhs.data_view();
   g.singularity_view() = rhs.singularity_view();
  }

 template<typename Descriptor, typename T>
  ENABLE_IF(arrays::is_scalar<T>) triqs_gf_view_assign_delegation( gf_view<Descriptor> & g, T const & x) {
   g._data_assigner_to_scalar(x, g.storage_is_array); //g.data_view() = x;
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
