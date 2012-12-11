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
#ifndef TRIQS_GF_GFBASECLASS_H
#define TRIQS_GF_GFBASECLASS_H
#include <triqs/utility/first_include.hpp>
#include "./tools.hpp"
#include <triqs/arrays/h5.hpp>
#include <vector>

namespace triqs { namespace gf { 

 template<typename Descriptor> class gf;         // the value class
 template<typename Descriptor> class gf_view;    // the view class

 // ---------------------- implementation --------------------------------

 // A class to store the values of the gf on a mesh when they are not array/matrix
 // Like a vector, but dressed ....
 
 template< typename T, bool IsView> class vector_storage { 
 public:
  typedef typename view_type_if_exists_else_type<T>::type  T_view_t;
  typedef typename non_view_type_if_exists_else_type<T>::type  T_non_view_t;
  typedef typename mpl::if_c<IsView, T_view_t, T>::type    T_t;
  typedef typename mpl::if_c<!IsView, T_view_t, T>::type   T_other_t;
  
  typedef void has_view_type_tag;  // Idiom : ValueView  
  typedef vector_storage<T,true>  view_type;
  typedef vector_storage<T,false> non_view_type;

  vector_storage(){}

  vector_storage( vector_storage const &x) : data( x.data) {}
  vector_storage( vector_storage &&x)      : data( std::move(x.data)) {}
  vector_storage( vector_storage<T, !IsView> const & V) { data.reserve(V.data.size()); for (auto & x : V.data ) data.push_back(x); } 

  vector_storage( std::vector<T_t>  const &V) : data(V){ } 
  vector_storage( std::vector<T_t>  &&V) : data( std::move(V)) { }

  vector_storage( std::vector<T_other_t>  const &V){ data.reserve(V.size()); for (auto & x : V) data.push_back(x); }
  vector_storage( std::vector<T_other_t>  &&V)     { data.reserve(V.size()); for (auto & x : V) data.push_back(x); } 

  template<typename V2, bool _is_view = IsView>
   ENABLE_IFC(_is_view) rebind (V2 const & v) { data.clear(); for (auto & x : v ) data.push_back(x); }

  vector_storage & operator = ( vector_storage             const & V) { data = V.data; return *this;}
  vector_storage & operator = ( vector_storage<T, !IsView> const & V) { data.clear(); for (auto & x : V.data ) data.push_back(x); return *this;  }

  template< bool _is_view = IsView>
   ENABLE_IFC(_is_view) operator = (typename T::value_type const & y) { for (auto & x : data ) x=y; }

  T_view_t operator()(tqa::ellipsis, size_t i) const { return data[i];} // ellipsis just here to make later code simpler below...
  T_view_t operator[](size_t i) const { return data[i];} 

  typedef typename std::vector<T_t>::const_iterator const_iterator;
  const_iterator begin() const { return data.begin();}
  const_iterator end() const { return data.end();}

  bool is_empty() const { return data.size()==0;}

  size_t size() const { return data.size(); }

  friend vector_storage operator +(vector_storage const & A, vector_storage const & B) {
   if (A.size() != B.size()) TRIQS_RUNTIME_ERROR << " Trying to add 2 gf of different size";
   std::vector<T_non_view_t> r; r.reserve(A.size());
   for (size_t i=0; i<A.size(); ++i) r.push_back(A.data[i] + B.data[i]); 
   return vector_storage(r); 
  }

 private :
  friend class vector_storage<T,!IsView>; 
  std::vector<T_t> data;

  //  BOOST Serialization
  friend class boost::serialization::access;
  template<class Archive>
   void serialize(Archive & ar, const unsigned int version) {
    ar & boost::serialization::make_nvp("data",data);
   }
 };

 // 
 template<typename T> struct _basic_scalar { typedef T type; };
 template<typename D> struct _basic_scalar<gf_view<D>> : _basic_scalar<typename D::target_t::value_type >{};

 // to metacompute the storage type

 template <typename Target, typename SO, bool UsingBigArray> struct gf_storage_type;
 template <typename Target, typename SO> struct gf_storage_type<Target,SO,true>  
 { typedef arrays::array < typename Target::value_type,Target::rank +1,SO> type;};
 template <typename Target, typename SO> struct gf_storage_type<Target,SO,false> 
 { typedef vector_storage <Target, false> type;};

 //------------------------------------------------------------------------
 /// A common implementation class. Idiom : ValueView
 template<typename Descriptor,bool IsView> class gf_impl : Descriptor::tag {
  public : 

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

   typedef arrays::Option::C storage_order;
   //typedef arrays::Option::Fortran storage_order;

   constexpr static bool using_big_array_storage = arrays::is_amv_value_class<target_t>::value; // let descriptor decide ?

   typedef typename gf_storage_type<target_t, storage_order,using_big_array_storage>::type     data_non_view_t;
   typedef typename data_non_view_t::view_type                                   data_view_t;
   typedef typename mpl::if_c<IsView, data_view_t, data_non_view_t>::type        data_t;

   typedef typename mpl::if_c<IsView, singularity_view_t, singularity_non_view_t>::type   singularity_t;

   mesh_t const & mesh() const                 { return _mesh;}
   domain_t const & domain() const             { return _mesh.domain();}
   data_view_t data_view()                     { return data;} 
   const data_view_t data_view() const         { return data;}
   singularity_view_t singularity_view()             { return singularity;} 
   const singularity_view_t singularity_view() const { return singularity;}

   symmetry_t const & symmetry() const { return _symmetry;}
   indices_t const & indices() const { return _indices;}

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
   gf_impl(gf_impl const & x)                    : _mesh(x._mesh), data(x.data), singularity(x.singularity), _symmetry(x._symmetry), _indices(x._indices){}
   gf_impl(gf_impl<Descriptor,!IsView> const & x): _mesh(x.mesh()), data(x.data_view()), singularity(x.singularity_view()), _symmetry(x.symmetry()), _indices(x.indices()){} 

   // from the data directly
   gf_impl(mesh_t const & m, data_view_t const & dat, singularity_view_t const & ad, symmetry_t const & s, indices_t const & ind ) : 
    _mesh(m), data(dat), singularity(ad),_symmetry(s), _indices(ind){}

   gf_impl(mesh_t const & m, data_t && dat, singularity_view_t const & ad, symmetry_t const & s, indices_t const & ind ) : 
    _mesh(m), data(std::move(dat)), singularity(ad),_symmetry(s), _indices(ind) {}

   void operator = (gf_impl const & rhs) = delete; // done in derived class.

   bool is_empty() const { return data.is_empty();}

  public:

   // First, a simple () returns a view, like for an array...
   view_type operator()() const { return *this;}
   
   /// Calls are (perfectly) forwarded to the Descriptor::operator(), except mesh_point_t and when 
   /// there is at least one lazy argument ...
   template<typename Arg0, typename... Args >    // match any argument list, picking out the first type : () is not permitted
    typename std::add_const<
    typename boost::lazy_disable_if<  // disable the template if one the following conditions it true 
    boost::mpl::or_< // starting condition [OR] 
    boost::is_base_of< typename std::remove_reference<Arg0>::type, mesh_point_t>,  // Arg0 is (a & or a &&) to a mesh_point_t 
    clef::one_is_lazy<Arg0, Args...>                          // One of Args is a lazy expression
     >,                                                       // end of OR 
    std::result_of<typename Descriptor::evaluator(mesh_t, data_t, singularity_t, Arg0, Args...)> // what is the result type of call
     >::type     // end of lazy_disable_if
    >::type // end of add_Const 
     operator() (Arg0&& arg0, Args&&... args) const { return _evaluator(_mesh, data, singularity, std::forward<Arg0>( arg0), std::forward<Args>(args)...); }

   // Interaction with the CLEF library : calling the gf with any clef expression as argument build a new clef expression
   template<typename Arg0, typename ...Args>
    typename boost::lazy_enable_if<    // enable the template if
    clef::one_is_lazy<Arg0, Args...>,  // One of Args is a lazy expression
    std::result_of<clef::lazy_call<view_type>(Arg0, Args...)>
     >::type     // end of lazy_enable_if 
     operator()(Arg0 arg0, Args... args) const { 
      static_assert(sizeof...(Args) == Descriptor::arity-1, "Incorrect number of variable in call");
      return clef::lazy_call<view_type>(*this)(arg0,args...);
     }

   /// A direct access to the grid point 
   target_view_t operator() (mesh_point_t const & x)       { return data(tqa::ellipsis(),x.m->index_to_linear(x.index));}

   /// A direct access to the grid point (const version)
   target_view_t operator() (mesh_point_t const & x) const { return data(tqa::ellipsis(),x.m->index_to_linear(x.index));}

   /// A direct access to the grid point 
   template<typename... Args> 
    target_view_t on_mesh (Args&&... args)       { return data(tqa::ellipsis(),_mesh.index_to_linear(mesh_index_t(std::forward<Args>(args)...)));}

   /// A direct access to the grid point (const version)
   template<typename... Args> 
    const target_view_t on_mesh (Args&&... args) const { return data(tqa::ellipsis(),_mesh.index_to_linear(mesh_index_t(std::forward<Args>(args)...)));}

  private:
   struct _on_mesh_wrapper {
    gf_impl const & f; _on_mesh_wrapper (gf_impl const & _f) : f(_f) {}
    template <typename... Args> target_view_t operator ()(Args && ... args) const { return f.on_mesh(std::forward<Args>(args)...);} 
   };
   _on_mesh_wrapper friend on_mesh(gf_impl const & f) { return f;}

  public:

   // Interaction with the CLEF library : auto assignment of the gf (gf(om_) = expression fills the functions by evaluation of expression)
   TRIQS_CLEF_HAS_AUTO_ASSIGN(); 
   template<typename F> friend void triqs_nvl_auto_assign (gf_impl & x, F f) { Descriptor::assign_from_expression(x._mesh,x.data, x.singularity,f);}

   /// [] Calls are (perfectly) forwarded to the Descriptor::operator[]
   //except when there is at least one lazy argument ...
   template<typename Arg >   
    typename boost::lazy_disable_if<  // disable the template if one the following conditions it true 
    clef::one_is_lazy<Arg>,                          // One of Args is a lazy expression
    std::result_of<typename Descriptor::bracket_evaluator(mesh_t, data_t &, singularity_t &, Arg)> // what is the result type of call
     >::type     // end of lazy_disable_if 
     operator[] (Arg&& arg) {return _bracket_evaluator(_mesh, this->data, singularity, std::forward<Arg>( arg));}

   /// [] Calls are (perfectly) forwarded to the Descriptor::operator[]
   //except when there is at least one lazy argument ...
   template<typename Arg >   
    typename boost::lazy_disable_if<  // disable the template if one the following conditions it true 
    clef::one_is_lazy<Arg>,                          // One of Args is a lazy expression
    std::result_of<typename Descriptor::bracket_evaluator(mesh_t, data_t &, singularity_t&, Arg)> // what is the result type of call
     >::type     // end of lazy_disable_if 
     operator[] (Arg&& arg) const {return _bracket_evaluator(_mesh, data, singularity, std::forward<Arg>( arg));}

   // Interaction with the CLEF library : calling the gf with any clef expression as argument build a new clef expression
   template<typename Arg>
    typename boost::lazy_enable_if<    // enable the template if
    clef::one_is_lazy<Arg>,  // One of Args is a lazy expression
    std::result_of<clef::lazy_call<view_type>(Arg)>
     >::type     // end of lazy_enable_if 
     operator[](Arg arg) const { return clef::lazy_call<view_type>(*this)(arg); }

   /// Write into HDF5
   friend void h5_write (tqa::h5::group_or_file fg, std::string subgroup_name, gf_impl const & g) {
    tqa::h5::group_or_file gr =  fg.create_group(subgroup_name);
    // Add the attribute
    h5_write(gr,"data",g.data);
    h5_write(gr,"singularity",g.singularity);
    h5_write(gr,"mesh",g._mesh);
    h5_write(gr,"symmetry",g._symmetry);
    h5_write(gr,"indices",g._indices);
   }

   /// Read from HDF5
   friend void h5_read  (tqa::h5::group_or_file fg, std::string subgroup_name, gf_impl & g){
    tqa::h5::group_or_file gr = fg.open_group(subgroup_name);
    // Check the attribute or throw 
    h5_read(gr,"data",g.data);
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

  gf(gf_view<Descriptor> const & g): B(g){} 

  template<typename GfType> gf(GfType const & x): B() { *this = x;} 

  gf(typename B::mesh_t const & m, 
    typename B::data_view_t const & dat, 
    typename B::singularity_view_t const & si, 
    typename B::symmetry_t const & s , 
    typename B::indices_t const & ind = typename B::indices_t () ) : 
   B(m,dat,si,s, ind) {}

  void operator = (gf const & rhs) {
   this->_mesh = rhs.mesh(); this->data = rhs.data; this->singularity = rhs.singularity;
   this->_symmetry = rhs._symmetry; this->_indices = rhs._indices;
  }

  template<typename RHS> void operator = (RHS const & rhs) { 
   this->_mesh = rhs.mesh(); this->data = rhs.data_view(); this->singularity = rhs.singularity_view();
   // There is a pb hrer : gf_proto has no symmetry .....
   //   this->_symmetry = rhs.symmetry(); this->_indices = rhs._indices();
  }
 };

 // ---------------------------------------------------------------------------------
 ///The View class of GF
 template<typename Descriptor> class gf_view : public gf_impl<Descriptor,true> {
  typedef gf_impl<Descriptor,true> B;

#ifdef TRIQS_ALLOW_EMPTY_VIEW
  public:
#endif
  gf_view ():B(){}

  public :

  void rebind( gf_view const &X) { 
   this->_mesh = X._mesh; this->_symmetry = X._symmetry; this->_indices = X._indices; 
   this->data.rebind(X.data); this->singularity.rebind(X.singularity); 
  }

  gf_view(gf_view const & g): B(g){}
  template<bool V> gf_view(gf_impl<Descriptor,V> const & g): B(g){}

  template<typename D, typename T> 
   gf_view (typename B::mesh_t const & m, 
     D const & dat,T const & t,typename B::symmetry_t const & s, 
     typename B::indices_t const & ind = typename B::indices_t () ) : B(m,dat,t,s,ind) {}

  void operator = (gf_view const & rhs)  { if (this->is_empty()) { rebind(rhs); return;} triqs_gf_view_assign_delegation(*this,rhs);}

  template<typename RHS> void operator = (RHS const & rhs) { triqs_gf_view_assign_delegation(*this,rhs);}
 };

 // delegate = so that I can overload it for specific RHS...
 template<typename Descriptor, typename RHS> 
  void triqs_gf_view_assign_delegation( gf_view<Descriptor> & g, RHS && rhs) { 
   if (!(g.mesh()  == rhs.mesh()))  TRIQS_RUNTIME_ERROR<<"Gf Assignment in View : incompatible mesh";
   //if (!(g.shape() == rhs.shape())) TRIQS_RUNTIME_ERROR<<"Gf Assignment in View : incompatible target shape";
   g.data_view() = rhs.data_view();
   g.singularity_view() = rhs.singularity_view();
  }

 template<typename Descriptor, typename T>
  ENABLE_IF(arrays::is_scalar<T>) triqs_gf_view_assign_delegation( gf_view<Descriptor> & g, T const & x) { g.data_view() = x; g.singularity_view() = x;}

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
#endif
