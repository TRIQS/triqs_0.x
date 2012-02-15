
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

#ifndef TRIQS_ARRAY_VECTOR_H
#define TRIQS_ARRAY_VECTOR_H
#include "indexmaps/cuboid/cuboid_map.hpp"
#include "indexmaps/cuboid/cuboid_slice.hpp"
#include "impl/indexmap_storage_pair.hpp"
#include "impl/providers.hpp"
#include "impl/option.hpp"

namespace triqs { namespace arrays {

 template <typename ValueType, typename Opt = Option::Default> class vector_view;
 template <typename ValueType, typename Opt = Option::Default> class vector;

 /** */
 template <typename ValueType, typename Opt >
  class vector_view : Tag::vector_view, Tag::vector_algebra_expression_terminal,
  public details::indexmap_storage_pair < typename R_Opt_2_IM<1,Opt>::type, storages::shared_block<ValueType>, Opt, Tag::vector_view >,
  public providers::compound_assign_ops<vector_view<ValueType,Opt> >
 {
  public :
   typedef details::indexmap_storage_pair < typename R_Opt_2_IM<1,Opt>::type, storages::shared_block<ValueType>, Opt, Tag::vector_view > BaseType;
   typedef vector_view<ValueType,Opt> view_type;
   typedef vector<ValueType,Opt> non_view_type;

   /// Build from an IndexMap and a storage 
   template<typename S> vector_view (indexmaps::cuboid_map<indexmaps::IndexOrder::C<1>, false > const & Ind,S const & Mem): BaseType(Ind, Mem) {}

   /// Build from anything that has an indexmap and a storage compatible with this class
   template<typename ISP>
    vector_view(const ISP & X): BaseType(X.indexmap(),X.storage()) {}

#ifdef TRIQS_ARRAYS_WITH_PYTHON_SUPPORT
   /**
    * Build from a numpy : only if TRIQS_ARRAYS_WITH_PYTHON_SUPPORT is defined
    */
   explicit vector_view (PyObject * X); // implemented in python/numpy_interface : include before use
#endif

   /// Copy construction
   vector_view(vector_view const & X): BaseType(X.indexmap(),X.storage()) {}

   /** Assignment.  The size of the array MUST match exactly.  */
   template<typename RHS> vector_view & operator=(const RHS & X) { assignment(*this,X); return *this; }

   vector_view & operator=(vector_view const & X) { assignment(*this,X); return *this; }//cf array_view class comment

   size_t size() const { return this->shape()[0];} 
   
 };

 template < class V, int R, class Opt > struct ViewFactory< V, R, Opt, Tag::vector_view> { typedef vector_view<V,Opt> type; };

 template <typename ValueType, typename Opt>
  class vector: Tag::vector,Tag::vector_algebra_expression_terminal,
  public vector_view<ValueType,Opt>::BaseType,
  public providers::compound_assign_ops<vector_view<ValueType,Opt> >
 {
  public :
   typedef typename vector_view<ValueType,Opt>::BaseType  BaseType;
   typedef typename BaseType::value_type value_type;
   typedef typename BaseType::storage_type storage_type;
   typedef typename BaseType::indexmap_type indexmap_type;
   typedef vector_view<ValueType,Opt> view_type;
   typedef vector<ValueType,Opt> non_view_type;

   /// Empty vector.
   vector():BaseType(indexmap_type(),storage_type()) {}

   ///
   vector(size_t dim):BaseType(indexmap_type(mini_vector<size_t,1>(dim))) {}

#ifdef TRIQS_ARRAYS_WITH_PYTHON_SUPPORT
   /**
    * Build from a numpy : only if TRIQS_ARRAYS_WITH_PYTHON_SUPPORT is defined
    */
   explicit vector (PyObject * X); // implemented in python/numpy_interface : include before use
#endif

   /** Makes a true (deep) copy of the data. */
   vector(const vector & X): BaseType(X.indexmap(),X.storage().clone()) {}

   /** 
    * Build a new vector from X.domain() and fill it with by evaluating X. X can be : 
    *  - another type of array, array_view, matrix,.... (any <IndexMap, Storage> pair)
    *  - a expression : e.g. array<int, IndexOrder::C<1> > A( B+ 2*C);
    */
   template <typename T> 
    vector(const T & X, typename boost::enable_if< is_array_assign_lhs<T> >::type *dummy =0):
     BaseType(indexmap_type(X.domain())) { assignment(*this,X); }

   /** 
    * Resizes the vector. NB : all references to the storage is invalidated.
    * Does not initialize the vector by default: to resize and init, do resize(IND).init()
    */
   vector & resize (size_t L) { BaseType::resize(typename BaseType::domain_type(mini_vector<size_t,1>(L))); return *this; }

   /// Assignement resizes the vector.  All references to the storage are therefore invalidated.
   vector & operator=(const vector & X) { BaseType::resize_and_clone_data(X); return *this; }

   /** 
    * Assignement resizes the vector.  All references to the storage are therefore invalidated.
    * NB : to avoid that, do make_view(A) = X instead of A = X
    */
   template<typename RHS> 
    vector & operator=(const RHS & X) { 
     static_assert(  is_array_assign_lhs<RHS>::value, "Assignment : RHS not supported");
     BaseType::resize(X.domain());
     assignment(*this,X);
     return *this; 
    }

   size_t size() const { return this->shape()[0];} 

 };//vector class
}}//namespace triqs::arrays

#include <boost/numeric/bindings/detail/adaptor.hpp>
namespace boost { namespace numeric { namespace bindings { namespace detail {

 template <typename ValueType, typename Opt, typename Id >
  struct adaptor<  triqs::arrays::vector_view<ValueType,Opt>, Id > {  
   typedef typename copy_const< Id, ValueType >::type value_type;
   typedef mpl::map<
    mpl::pair< tag::value_type, value_type >,
    mpl::pair< tag::entity, tag::vector >,
    mpl::pair< tag::size_type<1>, std::ptrdiff_t >,
    mpl::pair< tag::data_structure, tag::linear_array >,
    mpl::pair< tag::stride_type<1>,  std::ptrdiff_t > //tag::contiguous >
     > property_map;
   static std::ptrdiff_t size1( const Id& t ) { return t.indexmap().domain().number_of_elements(); }
   static value_type* begin_value( Id& t ) { return &(t(0)); }
   static value_type* end_value( Id& t ) { return &(t(t.indexmap().number_of_elements()-1)) + 1;}
   static std::ptrdiff_t stride1( const Id& id ) { return id.indexmap().strides()[0]; }
  };

 template <typename ValueType, typename Opt, typename Id >
  struct adaptor< triqs::arrays::vector<ValueType,Opt>, Id >: 
  adaptor<  triqs::arrays::vector_view<ValueType,Opt>, Id > {};

}}}} // namespace detail, namespace binding, namespace numeric, namespace boost


#include <boost/numeric/bindings/blas/level1/scal.hpp>
#include <boost/numeric/bindings/blas/level1/copy.hpp>
#include <boost/numeric/bindings/blas/level1/axpy.hpp>
#include <boost/numeric/bindings/blas/level1/swap.hpp>
namespace triqs { namespace arrays { namespace details { 

 // = for vector & vector_view.
 template<typename LHS, typename RHS> 
  struct assign_impl<LHS ,RHS, typename boost::enable_if<boost::mpl::and_<is_vector_or_view<LHS>,is_vector_or_view <RHS > > >::type > { 
   LHS & lhs; const RHS & rhs;
   assign_impl(LHS & lhs_, const RHS & rhs_): lhs(lhs_), rhs(rhs_) {}
   void invoke() { // std::cerr<<" vector view = with blas"<<std::endl;
    boost::numeric::bindings::blas::copy(rhs,lhs);} 
  };

 // *= and /= for scalar RHS, for vector_view.
 template<typename LHS, typename RHS, char OP> 
  struct comp_assign_impl<LHS,RHS,OP, 
  typename boost::enable_if<boost::mpl::and_<_is_MD<OP>, is_scalar_for<RHS,LHS >, is_vector_or_view<LHS> > >::type > { 
   LHS & lhs; const RHS & rhs;
   comp_assign_impl(LHS & lhs_, const RHS & rhs_): lhs(lhs_), rhs(rhs_) {}
   void invoke() {  // std::cerr<<" vector *= with blas"<<std::endl;
    typename LHS::value_type  a = (OP=='M' ? rhs : 1/rhs); boost::numeric::bindings::blas::scal(a,lhs);} 
  };

 // += and -= for vector with blas
 template<typename LHS, typename RHS,char OP> 
  struct comp_assign_impl<LHS,RHS,OP, 
  typename boost::enable_if<boost::mpl::and_<_is_AS<OP>, is_vector_or_view <RHS >, is_vector_or_view<LHS>  > >::type > { 
   LHS & lhs; const RHS & rhs;
   comp_assign_impl(LHS & lhs_, const RHS & rhs_): lhs(lhs_), rhs(rhs_) {}
   void invoke() { //  std::cerr<<" vector += with blas"<<std::endl;
    typename LHS::value_type a = (OP=='A' ? 1.0 : -1.0); boost::numeric::bindings::blas::axpy(a,rhs,lhs);} 
  };

}

// swapping 2 vector 
template <typename V, typename S1, typename S2>
void swap(vector_view <V,S1> x, vector_view<V,S2> y) { boost::numeric::bindings::blas::swap(x,y);} 

template <typename V, typename S1, typename S2>
void swap(vector <V,S1> & x, vector<V,S2>  & y) { boost::numeric::bindings::blas::swap(x,y);} 

}}
#endif

