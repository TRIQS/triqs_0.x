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
#include "indexmaps/cuboid/cuboid_map.hpp"
#include "indexmaps/cuboid/cuboid_slice.hpp"
#include "impl/indexmap_storage_pair.hpp"
#include "impl/assignment.hpp"
#include "vector.hpp"

namespace triqs { namespace arrays {

 template <typename ValueType, typename Opt = Option::Default> class matrix_view;
 template <typename ValueType, typename Opt = Option::Default> class matrix;

 namespace details { // ********************************************************
  using namespace indexmaps; 

  template<typename OrderTag> struct order_as_char;
  template<> struct order_as_char<Tag::C> { static const char val='C';};
  template<> struct order_as_char<Tag::Fortran> { static const char val='F';};
  // add the  permutations ?

  // flip the order, for transposition.
  template<char C> struct flip_order;
  template<> struct flip_order<'C'> { typedef Tag::Fortran type;};
  template<> struct flip_order<'F'> { typedef Tag::C type;};

  template<typename IO> struct tag_from_index_order;
  template<int D> struct tag_from_index_order<IndexOrder::C<D> > { typedef Tag::C type; };
  template<int D> struct tag_from_index_order<IndexOrder::Fortran<D> > { typedef Tag::Fortran type; };

 }// details ********************************************************

 // this traits is used by indexmap_storage_pair, when slicing to find the correct view type.
 template < class V, int R, class Opt > struct ViewFactory< V, R, Opt, Tag::matrix_view> { 
  typedef typename boost::mpl::if_c<R == 1, vector_view<V,Opt >, matrix_view<V,Opt > >::type type;
 };

 #define _IMPL_MATRIX_COMMON \
  size_t dim0() const { return this->shape()[0];}\
  size_t dim1() const { return this->shape()[1];}\
  bool is_square() const { return dim0() == dim1();}\
  \
  static const char order = details::order_as_char<typename Opt::IndexOrderTag>::val;\
  typedef matrix_view<ValueType, Option::options<typename details::flip_order<order>::type, typename Opt::StorageTag > > transpose_return_type;\
  transpose_return_type transpose() const {\
   typedef typename transpose_return_type::indexmap_type im_type; typedef Permutations::permutation<1,0> P;\
   return transpose_return_type( im_type(eval_permutation<P>(this->indexmap().lengths()), eval_permutation<P>(this->indexmap().strides()),\
      this->indexmap().start_shift()), this->storage());\
  }

#define IMPL_TYPE details::indexmap_storage_pair < typename R_Opt_2_IM<2,Opt>::type, storages::shared_block<ValueType>, Opt, Tag::matrix_view > 

 template <typename ValueType, typename Opt >
  class matrix_view : Tag::matrix_view,  TRIQS_MODEL_CONCEPT(ImmutableMatrix), public IMPL_TYPE {
   public :
   typedef matrix_view<ValueType,Opt> view_type;
   typedef matrix<ValueType,Opt>      non_view_type;
   typedef void has_view_type_tag;
   typedef Opt opt_type;

   /// Build from an IndexMap and a storage 
   template<typename S> matrix_view (typename IMPL_TYPE::indexmap_type const & Ind,S const & Mem): IMPL_TYPE(Ind, Mem) {}

   /// Build from anything that has an indexmap and a storage compatible with this class
   template<typename ISP> matrix_view(const ISP & X): IMPL_TYPE(X.indexmap(),X.storage()) {}

#ifdef TRIQS_ARRAYS_WITH_PYTHON_SUPPORT
   /// Build from a numpy.array : throws if X is not a numpy.array 
   explicit matrix_view (PyObject * X): IMPL_TYPE(X, false, "matrix_view "){}

   /// Build from a numpy.array : throws if X is not a numpy.array 
   explicit matrix_view (boost::python::object X): IMPL_TYPE(X.ptr(), false, "matrix_view "){}
#endif

   /// Copy construction
   matrix_view( matrix_view const & X): IMPL_TYPE(X.indexmap(),X.storage()) {}

   /** Assignement.  The size of the array MUST match exactly.  */
   template<typename RHS> matrix_view & operator=(const RHS & X) {triqs_arrays_assign_delegation(*this,X); return *this; }

   matrix_view & operator=(matrix_view const & X) { triqs_arrays_assign_delegation(*this,X); return *this; }//cf array_view class comment

   TRIQS_DEFINE_COMPOUND_OPERATORS(matrix_view); 
   _IMPL_MATRIX_COMMON;
   };

 //--------------------------------------------------

 template <typename ValueType, typename Opt >
  class matrix: Tag::matrix,  TRIQS_MODEL_CONCEPT(ImmutableMatrix), public IMPL_TYPE {
   public :
   typedef typename IMPL_TYPE::value_type value_type;
   typedef typename IMPL_TYPE::storage_type storage_type;
   typedef typename IMPL_TYPE::indexmap_type indexmap_type;
   typedef matrix_view<ValueType,Opt> view_type;
   typedef matrix<ValueType,Opt> non_view_type;
   typedef void has_view_type_tag;
   typedef Opt opt_type;

   /// Empty matrix.
   matrix():IMPL_TYPE(indexmap_type(),storage_type()) {}

   ///
   matrix(size_t dim1, size_t dim2) : IMPL_TYPE(indexmap_type(mini_vector<size_t,2>(dim1,dim2))) {}
   
   ///
   matrix(mini_vector<size_t,2> const & sha) : IMPL_TYPE(indexmap_type(sha)) {}

   /** Makes a true (deep) copy of the data. */
   matrix(const matrix & X): IMPL_TYPE(X.indexmap(),X.storage().clone()) {}

   /// Build a new matrix from X.domain() and fill it with by evaluating X. X can be : 
   template <typename T> 
    matrix(const T & X, typename boost::enable_if< ImmutableArray<T> >::type *dummy =0):
     IMPL_TYPE(indexmap_type(X.domain())) { triqs_arrays_assign_delegation(*this,X); }

#ifdef TRIQS_ARRAYS_WITH_PYTHON_SUPPORT
   ///Build from a numpy.array X (or any object from which numpy can make a numpy.array). Makes a copy.
   explicit matrix (PyObject * X): IMPL_TYPE(X, true, "matrix "){}

   ///Build from a numpy.array X (or any object from which numpy can make a numpy.array). Makes a copy.
   explicit matrix (boost::python::object X): IMPL_TYPE(X.ptr(), true, "matrix "){}
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
   matrix & resize (const indexmaps::cuboid_domain<IMPL_TYPE::rank> & l) { IMPL_TYPE::resize(l); return *this; }

   /// Assignement resizes the matrix.  All references to the storage are therefore invalidated.
   matrix & operator=(const matrix & X) { IMPL_TYPE::resize_and_clone_data(X); return *this; }

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
 namespace details { 

  template <typename T, int R> 
   bool kronecker(mini_vector<T,R> const & key) { return ( (R==2) && (key[0]==key[1]));} 

  // assignment for scalar RHS // write a specific one if it is not a view : plain loop
  // beware : for matrix, assign to a scalar will make the matrix scalar, as it should

  template <typename V, typename RHS> struct __looper { 
   RHS rhs;
   __looper(const RHS & rhs_): rhs(rhs_) {}
   template <typename KeyType> void operator()(V & p, KeyType const & key) const { p = (kronecker(key) ? rhs : RHS() ); }
  };

  template<typename RHS, typename V, typename Opt> 
   typename boost::enable_if<is_scalar_for<RHS,matrix_view<V,Opt> > >::type
   triqs_arrays_assign_delegation (matrix<V,Opt> & lhs, RHS const & rhs) { indexmaps::foreach( __looper<V,RHS>(rhs),lhs); }

  template<typename RHS, typename V, typename Opt> 
   typename boost::enable_if<is_scalar_for<RHS,matrix_view<V,Opt> > >::type
   triqs_arrays_assign_delegation (matrix_view<V,Opt> & lhs, RHS const & rhs) { indexmaps::foreach( __looper<V,RHS>(rhs),lhs); }

  /* 
     template <typename V, typename Opt, typename RHS> 
     struct assign_impl<matrix_view<V,Opt>,RHS,typename boost::enable_if<is_scalar_for<RHS,matrix_view<V,Opt> > >::type > { 
     typedef matrix_view<V,Opt> LHS;
     typedef typename LHS::indexmap_type::domain_type::index_value_type index_value_type;
     LHS & lhs; const RHS & rhs;
     assign_impl(LHS & lhs_, const RHS & rhs_): lhs(lhs_), rhs(rhs_) {}
     void operator()(V & p, index_value_type const & key) const { p = (kronecker(key) ? rhs : RHS() ); }
     void invoke() { foreach(*this,lhs); } // if contiguous : plain loop else foreach...
     };

     template <typename V, typename Opt, typename RHS> 
     struct assign_impl<matrix<V,Opt>,RHS,typename boost::enable_if<is_scalar_for<RHS,matrix_view<V,Opt> > >::type > : 
     assign_impl<matrix_view<V,Opt>,RHS > {
     typedef matrix<V,Opt> LHS;
     assign_impl(LHS & lhs_, const RHS & rhs_): assign_impl<matrix_view<V,Opt>,RHS > (lhs_,rhs_){}
     }; 
     */

 } //details
}}//namespace triqs::arrays

#include <boost/numeric/bindings/detail/adaptor.hpp>
#include <boost/numeric/bindings/detail/if_row_major.hpp>
// blas lapack binder
namespace boost { namespace numeric { namespace bindings { namespace detail {

 namespace numerical_array_detail { 
  template <char C> struct order_trait;
  template <> struct order_trait<'F'> { typedef  tag::column_major data_order;};
  template <> struct order_trait<'C'> { typedef  tag::row_major data_order;};
 }

 template <typename ValueType, typename Opt, typename Id >
  struct adaptor<  triqs::arrays::matrix_view <ValueType,Opt>, Id > {
   typedef typename copy_const< Id, ValueType >::type value_type;
   typedef typename  numerical_array_detail::order_trait< triqs::arrays::matrix_view <ValueType,Opt>::order >::data_order  data_order;

   typedef mpl::map<
    mpl::pair< tag::value_type, value_type >,
    mpl::pair< tag::entity, tag::matrix >,
    mpl::pair< tag::size_type<1>, std::ptrdiff_t >,
    mpl::pair< tag::size_type<2>, std::ptrdiff_t >,
    mpl::pair< tag::data_structure, tag::linear_array >,
    mpl::pair< tag::data_order, data_order  >,
    mpl::pair< tag::stride_type<1>,
    typename if_row_major< data_order, std::ptrdiff_t, tag::contiguous >::type >,
    mpl::pair< tag::stride_type<2>,
    typename if_row_major< data_order, tag::contiguous, std::ptrdiff_t >::type >
     > property_map;

   static std::ptrdiff_t size1( const Id& id ) { return id.indexmap().lengths()[0]; }
   static std::ptrdiff_t size2( const Id& id ) { return id.indexmap().lengths()[1];}
   static value_type* begin_value( Id& t ) { return &(t(0,0)); }
   static value_type* end_value( Id& t ) { return &t(t.dim0()-1,t.dim1()-1) + 1; } // can't find the doc : if +1 needed ???
   static std::ptrdiff_t stride1( const Id& id ) { return id.indexmap().strides()[0]; }
   static std::ptrdiff_t stride2( const Id& id ) { return id.indexmap().strides()[1]; }
  };

 template <typename ValueType, typename Opt, typename Id >
  struct adaptor< triqs::arrays::matrix<ValueType,Opt>, Id >: adaptor<  triqs::arrays::matrix_view<ValueType,Opt>, Id > {};

}}}} // namespace detail, namespace binding, namespace numeric, namespace boost
#endif

