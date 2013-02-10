/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011 by O. Parcollet
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
#ifndef TRIQS_ARRAYS_MATRIX_H
#define TRIQS_ARRAYS_MATRIX_H
#include "indexmaps/cuboid/map.hpp"
#include "indexmaps/cuboid/slice.hpp"
#include "impl/indexmap_storage_pair.hpp"
#include "impl/assignment.hpp"
#include "vector.hpp"

namespace triqs { namespace arrays {

 template <typename ValueType, ull_t Opt=0> class matrix_view;
 template <typename ValueType, ull_t Opt=0> class matrix;
 
 // ---------------------- matrix --------------------------------
 //
#define _IMPL_MATRIX_COMMON \
 size_t dim0() const { return this->shape()[0];}\
 size_t dim1() const { return this->shape()[1];}\
 bool is_square() const { return dim0() == dim1();}\
 \
 view_type transpose() const {\
  typename indexmap_type::lengths_type l; l[0] = this->indexmap().lengths()[1];l[1] = this->indexmap().lengths()[0];\
  typename indexmap_type::strides_type s; s[0] = this->indexmap().strides()[1];s[1] = this->indexmap().strides()[0];\
  return view_type( indexmap_type(l,s, this->indexmap().start_shift()), this->storage());\
 }

#define IMPL_TYPE indexmap_storage_pair < indexmaps::cuboid::map<2,Opt>, storages::shared_block<ValueType>, Opt, Tag::matrix_view > 

 template <typename ValueType, ull_t Opt >
  class matrix_view : Tag::matrix_view,  TRIQS_MODEL_CONCEPT(ImmutableMatrix), public IMPL_TYPE {
   public :
    typedef matrix_view<ValueType,Opt> view_type;
    typedef matrix<ValueType,Opt>      non_view_type;
    typedef void has_view_type_tag;

    typedef typename IMPL_TYPE::indexmap_type indexmap_type;
    typedef typename IMPL_TYPE::storage_type storage_type;

    /// Build from an IndexMap and a storage 
    template<typename S> matrix_view (typename IMPL_TYPE::indexmap_type const & Ind,S const & Mem): IMPL_TYPE(Ind, Mem) {}

    /// Build from anything that has an indexmap and a storage compatible with this class
    template<typename ISP> matrix_view(const ISP & X): IMPL_TYPE(X.indexmap(),X.storage()) {}

#ifdef TRIQS_WITH_PYTHON_SUPPORT
    /// Build from a numpy.array : throws if X is not a numpy.array 
    explicit matrix_view (PyObject * X): IMPL_TYPE(X, false, "matrix_view "){}
#endif

    /// Copy construction
    matrix_view( matrix_view const & X): IMPL_TYPE(X.indexmap(),X.storage()) {}

    /// Rebind the view
    void rebind (matrix_view const & X) { this->indexmap_ = X.indexmap_; this->storage_ = X.storage_;}

    /** Assignement.  The size of the array MUST match exactly.  */
    template<typename RHS> matrix_view & operator=(const RHS & X) {triqs_arrays_assign_delegation(*this,X); return *this; }

    matrix_view & operator=(matrix_view const & X) { 
     if (this->is_empty()) rebind(X); else triqs_arrays_assign_delegation(*this,X); return *this; }//cf array_view class comment

    TRIQS_DEFINE_COMPOUND_OPERATORS(matrix_view); 
    _IMPL_MATRIX_COMMON;
  };

 //---------------------------------------------------------------------
 // this traits is used by indexmap_storage_pair, when slicing to find the correct view type.
 template < class V, int R, ull_t OptionFlags > struct ViewFactory< V, R, OptionFlags, Tag::matrix_view> { 
  typedef typename std::conditional <R == 1, vector_view<V,OptionFlags >, matrix_view<V,OptionFlags > >::type type;
 };

 // ---------------------- matrix --------------------------------

 template <typename ValueType, ull_t Opt >
  class matrix: Tag::matrix,  TRIQS_MODEL_CONCEPT(ImmutableMatrix), public IMPL_TYPE {
   public :
    typedef typename IMPL_TYPE::value_type value_type;
    typedef typename IMPL_TYPE::storage_type storage_type;
    typedef typename IMPL_TYPE::indexmap_type indexmap_type;
    typedef matrix_view<ValueType,Opt> view_type;
    typedef matrix<ValueType,Opt> non_view_type;
    typedef void has_view_type_tag;

    /// Empty matrix.
   // ambigous
   // matrix(char ml='C'):  IMPL_TYPE(indexmap_type(memory_layout<2>(ml))) {}

    /// Empty matrix.
    matrix(memory_layout<2> ml = memory_layout<2>() ):  IMPL_TYPE(indexmap_type(ml)) {}

    /// Move
    explicit matrix(matrix && X) { this->swap_me(X); } 

    ///
    matrix(size_t dim1, size_t dim2, memory_layout<2> ml = memory_layout<2>() ) : IMPL_TYPE(indexmap_type(mini_vector<size_t,2>(dim1,dim2),ml)) {}

    ///
    matrix(mini_vector<size_t,2> const & sha, memory_layout<2> ml = memory_layout<2>()) : IMPL_TYPE(indexmap_type(sha,ml)) {}

    /** Makes a true (deep) copy of the data. */
    matrix(const matrix & X): IMPL_TYPE(X.indexmap(),X.storage().clone()) {}

    /// Build a new matrix from X.domain() and fill it with by evaluating X. X can be : 
    template <typename T> 
     matrix(const T & X, TYPE_ENABLE_IF(memory_layout<2>, ImmutableArray<T>) ml = memory_layout<2>()):
     //matrix(const T & X, typename boost::enable_if< ImmutableArray<T> >::type *dummy =0):
      IMPL_TYPE(indexmap_type(X.domain(),ml)) { triqs_arrays_assign_delegation(*this,X); }

#ifdef TRIQS_WITH_PYTHON_SUPPORT
    ///Build from a numpy.array X (or any object from which numpy can make a numpy.array). Makes a copy.
    explicit matrix (PyObject * X): IMPL_TYPE(X, true, "matrix "){}
#endif

    /** 
     * Resizes the matrix. NB : all references to the storage is invalidated.
     * Does not initialize the matrix by default
     */
    matrix & resize (size_t n1, size_t n2) { IMPL_TYPE::resize(typename IMPL_TYPE::domain_type (mini_vector<size_t,2>(n1,n2))); return *this;}

    /** 
     * Resizes the matrix. NB : all references to the storage is invalidated.
     * Does not initialize the matrix by default
     */
    matrix & resize (const indexmaps::cuboid::domain<IMPL_TYPE::rank> & l) { IMPL_TYPE::resize(l); return *this; }

    /// Assignement resizes the matrix.  All references to the storage are therefore invalidated.
    matrix & operator=(const matrix & X) { IMPL_TYPE::resize_and_clone_data(X); return *this; }

    /// Move assignment
    matrix & operator=(matrix && X) { swap(*this, X); return *this;}

    /** 
     * Assignement resizes the matrix.  All references to the storage are therefore invalidated.
     * NB : to avoid that, do make_view(A) = X instead of A = X
     */
    template<typename RHS> 
     matrix & operator=(const RHS & X) { 
      static_assert(  ImmutableArray<RHS>::value, "Assignment : RHS not supported");
      IMPL_TYPE::resize(X.domain());
      triqs_arrays_assign_delegation(*this,X);
      return *this; 
     }

    TRIQS_DEFINE_COMPOUND_OPERATORS(matrix);
    _IMPL_MATRIX_COMMON;
  };//matrix class
#undef _IMPL_MATRIX_COMMON
#undef IMPL_TYPE

 template <typename T, int R> 
  bool kronecker(mini_vector<T,R> const & key) { return ( (R==2) && (key[0]==key[1]));} 

 // assignment for scalar RHS // write a specific one if it is not a view : plain loop
 // beware : for matrix, assign to a scalar will make the matrix scalar, as it should

 template <typename V, typename RHS> struct __looper { 
  RHS rhs;
  __looper(const RHS & rhs_): rhs(rhs_) {}
  template <typename KeyType> void operator()(V & p, KeyType const & key) const { p = (kronecker(key) ? rhs : RHS() ); }
 };

 template<typename RHS, typename V, ull_t Opt> 
  typename boost::enable_if<is_scalar_for<RHS,matrix_view<V,Opt> > >::type
  triqs_arrays_assign_delegation (matrix<V,Opt> & lhs, RHS const & rhs) { indexmaps::foreach( __looper<V,RHS>(rhs),lhs); }

 template<typename RHS, typename V, ull_t Opt> 
  typename boost::enable_if<is_scalar_for<RHS,matrix_view<V,Opt> > >::type
  triqs_arrays_assign_delegation (matrix_view<V,Opt> & lhs, RHS const & rhs) { indexmaps::foreach( __looper<V,RHS>(rhs),lhs); }

 // += for scalar RHS // write a specific one if it is not a view : plain loop
 // beware : for matrix, assign to a scalar will make the matrix scalar, as it should
 template <typename V, typename RHS> struct __looper_add {
  RHS rhs;
  __looper_add(const RHS & rhs_): rhs(rhs_) {}
  template <typename KeyType> void operator()(V & p, KeyType const & key) const { p += (kronecker(key) ? rhs : RHS() ); }
 };

 template<typename RHS, typename V, ull_t Opt>
  typename boost::enable_if<is_scalar_for<RHS,matrix_view<V,Opt> > >::type
  triqs_arrays_compound_assign_delegation (matrix<V,Opt> & lhs, RHS const & rhs, mpl::char_<'A'>) { indexmaps::foreach( __looper_add<V,RHS>(rhs),lhs); }

 template<typename RHS, typename V, ull_t Opt>
  typename boost::enable_if<is_scalar_for<RHS,matrix_view<V,Opt> > >::type
  triqs_arrays_compound_assign_delegation (matrix_view<V,Opt> & lhs, RHS const & rhs, mpl::char_<'A'>) { indexmaps::foreach( __looper_add<V,RHS>(rhs),lhs); }

 template<typename RHS, typename V, ull_t Opt>
  typename boost::enable_if<is_scalar_for<RHS,matrix_view<V,Opt> > >::type
  triqs_arrays_compound_assign_delegation (matrix<V,Opt> & lhs, RHS const & rhs, mpl::char_<'S'>) { lhs += (-rhs); }

 template<typename RHS, typename V, ull_t Opt>
  typename boost::enable_if<is_scalar_for<RHS,matrix_view<V,Opt> > >::type
  triqs_arrays_compound_assign_delegation (matrix_view<V,Opt> & lhs, RHS const & rhs, mpl::char_<'S'>) {lhs += (-rhs); }

 template<typename RHS, typename V, typename Opt>
  typename std::enable_if<!is_scalar_for<RHS,matrix_view<V,Opt>>::value>::type
  triqs_arrays_compound_assign_delegation (matrix<V,Opt> & lhs, RHS const & rhs, mpl::char_<'M'>) = delete;

 template<typename RHS, typename V, typename Opt>
  typename std::enable_if<!is_scalar_for<RHS,matrix_view<V,Opt>>::value>::type
  triqs_arrays_compound_assign_delegation (matrix_view<V,Opt> & lhs, RHS const & rhs, mpl::char_<'M'>) = delete;

 template<typename RHS, typename V, typename Opt>
  typename std::enable_if<!is_scalar_for<RHS,matrix_view<V,Opt>>::value>::type
  triqs_arrays_compound_assign_delegation (matrix<V,Opt> & lhs, RHS const & rhs, mpl::char_<'D'>) = delete;

 template<typename RHS, typename V, typename Opt>
  typename std::enable_if<!is_scalar_for<RHS,matrix_view<V,Opt>>::value>::type
  triqs_arrays_compound_assign_delegation (matrix_view<V,Opt> & lhs, RHS const & rhs, mpl::char_<'D'>) = delete;


}}//namespace triqs::arrays

#include <boost/numeric/bindings/detail/adaptor.hpp>
#include <boost/numeric/bindings/detail/if_row_major.hpp>
// blas lapack binder
namespace boost { namespace numeric { namespace bindings { namespace detail {

 /*namespace numerical_array_detail { 
  template <char C> struct order_trait;
  template <> struct order_trait<'F'> { typedef  tag::column_major data_order;};
  template <> struct order_trait<'C'> { typedef  tag::row_major data_order;};
 }
*/
 template <typename ValueType, triqs::ull_t Opt, typename Id >
  struct adaptor<  triqs::arrays::matrix_view <ValueType,Opt>, Id > {
   typedef typename copy_const< Id, ValueType >::type value_type;
   //typedef typename  numerical_array_detail::order_trait< triqs::arrays::matrix_view <ValueType,Opt>::order >::data_order  data_order;
   
   // PB PB PB  : to be fixed : always row major !!!
   typedef tag::column_major data_order;

   typedef mpl::map<
    mpl::pair< tag::value_type, value_type >,
    mpl::pair< tag::entity, tag::matrix >,
    mpl::pair< tag::size_type<1>, std::ptrdiff_t >,
    mpl::pair< tag::size_type<2>, std::ptrdiff_t >,
    mpl::pair< tag::data_structure, tag::linear_array >,
    mpl::pair< tag::data_order, data_order  >,
    mpl::pair< tag::stride_type<1>,std::ptrdiff_t >,
    mpl::pair< tag::stride_type<2>,std::ptrdiff_t >
    //mpl::pair< tag::stride_type<1>,
    //typename if_row_major< data_order, std::ptrdiff_t, tag::contiguous >::type >,
    //mpl::pair< tag::stride_type<2>,
    //typename if_row_major< data_order, tag::contiguous, std::ptrdiff_t >::type >
     > property_map;

   //static std::ptrdiff_t size1( const Id& id ) { return id.indexmap().lengths()[0]; }
   //static std::ptrdiff_t size2( const Id& id ) { return id.indexmap().lengths()[1];}
   //static std::ptrdiff_t stride1( const Id& id ) { return id.indexmap().strides()[0]; }
   //static std::ptrdiff_t stride2( const Id& id ) { return id.indexmap().strides()[1]; }

   static bool is_f(Id const & id) { return  id.indexmap().memory_layout_is_fortran();}
   static std::ptrdiff_t size1( const Id& id ) { return id.indexmap().lengths()[is_f(id) ? 0 :1]; }
   static std::ptrdiff_t size2( const Id& id ) { return id.indexmap().lengths()[is_f(id) ? 1: 0];}
   static value_type* begin_value( Id& t ) { return &(t(0,0)); }
   static value_type* end_value( Id& t ) { return &t(t.dim0()-1,t.dim1()-1) + 1; } // can't find the doc : if +1 needed ???
   static std::ptrdiff_t stride1( const Id& id ) { return id.indexmap().strides()[is_f(id) ? 0 :1]; }
   static std::ptrdiff_t stride2( const Id& id ) { return id.indexmap().strides()[is_f(id) ? 1 :0]; }
  };

 template <typename ValueType, triqs::ull_t Opt, typename Id >
  struct adaptor< triqs::arrays::matrix<ValueType,Opt>, Id >: adaptor<  triqs::arrays::matrix_view<ValueType,Opt>, Id > {};

}}}} // namespace detail, namespace binding, namespace numeric, namespace boost

// The std::swap is WRONG for a view because of the copy/move semantics of view.
// Use swap instead (the correct one, found by ADL).
namespace std { 
 template <typename V, triqs::ull_t Opt >
  void swap( triqs::arrays::matrix_view<V,Opt> & a , triqs::arrays::matrix_view<V,Opt> & b)= delete;
}

#endif

