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
#ifndef TRIQS_ARRAYS_ARRAY_H
#define TRIQS_ARRAYS_ARRAY_H
#include "indexmaps/cuboid/cuboid_map.hpp"
#include "indexmaps/cuboid/cuboid_slice.hpp"
#include "impl/indexmap_storage_pair.hpp"
#include "impl/assignment.hpp"
#include "impl/option.hpp"
namespace triqs { namespace arrays {

 template <typename ValueType, int Rank, typename Opt= Option::Default > class array_view;
 template <typename ValueType, int Rank, typename Opt= Option::Default > class array;

 // ---------------------- implementation --------------------------------

#define IMPL_TYPE indexmap_storage_pair<typename R_Opt_2_IM<Rank,Opt>::type, storages::shared_block<ValueType>, Opt, Tag::array_view > 

 template <typename ValueType, int Rank, typename Opt>
  class array_view : Tag::array_view, TRIQS_MODEL_CONCEPT(ImmutableCuboidArray), public IMPL_TYPE {   
   static_assert( Rank>0, " Rank must be >0");
   public:   
   typedef typename IMPL_TYPE::indexmap_type indexmap_type;
   typedef typename IMPL_TYPE::storage_type storage_type;
   typedef array_view<ValueType,Rank,Opt> view_type;
   typedef array<ValueType,Rank,Opt> non_view_type;
   typedef void has_view_type_tag;

   /// Build from an IndexMap and a storage 
   template<typename S> array_view (indexmap_type const & Ind,S const & Mem): IMPL_TYPE(Ind, Mem) {}

   /// Copy constructor
   array_view(array_view const & X): IMPL_TYPE(X.indexmap(),X.storage()) {}

   /// Build from anything that has an indexmap and a storage compatible with this class
   template<typename ISP> array_view(const ISP & X): IMPL_TYPE(X.indexmap(),X.storage()) {}

#ifdef TRIQS_WITH_PYTHON_SUPPORT
   /// Build from a numpy.array : throws if X is not a numpy.array 
   explicit array_view (PyObject * X): IMPL_TYPE(X, false, "array_view "){}
#endif

#ifndef TRIQS_ALLOW_EMPTY_VIEW
   private:
#endif
   array_view (){}

   public:

   /// Rebind the view
   void rebind (array_view const & X) { this->indexmap_ = X.indexmap_; this->storage_ = X.storage_;}

   /// Assignment. The size of the array MUST match exactly, except in the empty case 
   template<typename RHS> array_view & operator=(RHS const & X) { triqs_arrays_assign_delegation(*this,X); return *this; }

   ///
   array_view & operator=(array_view const & X) {  
    if (this->is_empty()) rebind(X); else triqs_arrays_assign_delegation(*this,X); return *this; } //without this, the standard = is synthetized...

   TRIQS_DEFINE_COMPOUND_OPERATORS(array_view);
  };

 template < class V, int R, class Opt > struct ViewFactory< V, R, Opt, Tag::array_view > { typedef array_view<V,R,Opt> type; };

 template <typename ValueType, int Rank, typename Opt>
  class array: Tag::array,  TRIQS_MODEL_CONCEPT(ImmutableCuboidArray), public IMPL_TYPE { 
   public:
    typedef typename IMPL_TYPE::value_type value_type;
    typedef typename IMPL_TYPE::storage_type storage_type;
    typedef typename IMPL_TYPE::indexmap_type indexmap_type;
    typedef array_view<ValueType,Rank,Opt> view_type;
    typedef array<ValueType,Rank,Opt> non_view_type;
    typedef void has_view_type_tag;

    /// Empty array.
    array() {} 

    /// From a domain
    explicit array( typename indexmap_type::domain_type const & dom):IMPL_TYPE(indexmap_type(dom)){}

#ifdef TRIQS_DOXYGEN
    /// Construction from the dimensions. NB : the number of parameters must be exactly rank (checked at compile time). 
    array (size_t I_1, .... , size_t I_rank);
#else
#define IMPL(z, NN, unused)                                \
    explicit array (BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(NN), size_t I_)): \
    IMPL_TYPE(indexmap_type(mini_vector<size_t,BOOST_PP_INC(NN)>(BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(NN), I_)))) {\
     static_assert(IMPL_TYPE::rank-1==NN,"array : incorrect number of variables in constructor");}
    BOOST_PP_REPEAT(ARRAY_NRANK_MAX , IMPL, nil)
#undef IMPL
#endif

     // Makes a true (deep) copy of the data. 
     array(const array & X): IMPL_TYPE(X.indexmap(),X.storage().clone()) {}

     // Move
     explicit array(array && X) { this->swap_me(X); } 

    /** 
     * Build a new array from X.domain() and fill it with by evaluating X. X can be : 
     *  - another type of array, array_view, matrix,.... (any <IndexMap, Storage> pair)
     *  - a expression : e.g. array<int> A = B+ 2*C;
     */
    template <typename T> 
     array(const T & X, typename boost::enable_if< ImmutableArray<T> >::type *dummy =0):
      IMPL_TYPE(indexmap_type(X.domain())) { triqs_arrays_assign_delegation(*this,X); }

#ifdef TRIQS_WITH_PYTHON_SUPPORT
    ///Build from a numpy.array X (or any object from which numpy can make a numpy.array). Makes a copy.
    explicit array (PyObject * X): IMPL_TYPE(X, true, "array "){}
#endif

    /** 
     * Resizes the array. NB : all references to the storage is invalidated.
     * Does not initialize the array by default: to resize and init, do resize(IND).init()
     */
    array & resize (const indexmaps::cuboid_domain<IMPL_TYPE::rank> & l) { IMPL_TYPE::resize(l); return *this; }

    /// Assignement resizes the array.  All references to the storage are therefore invalidated.
    array & operator=(const array & X) { IMPL_TYPE::resize_and_clone_data(X); return *this; }

    /// Move assignment
    array & operator=(array && X) { swap(*this, X); return *this;}

    /** 
     * Assignement resizes the array (if necessary).
     * All references to the storage are therefore invalidated.
     * NB : to avoid that, do make_view(A) = X instead of A = X
     */
    template<typename RHS> 
     array & operator=(const RHS & X) { 
      static_assert(ImmutableArray<RHS>::value, "Assignment : RHS not supported");
      IMPL_TYPE::resize(X.domain());
      triqs_arrays_assign_delegation(*this,X);
      return *this; 
     }

    TRIQS_DEFINE_COMPOUND_OPERATORS(array);

  };//array class
}}//namespace triqs::arrays

#undef IMPL_TYPE

// The std::swap is WRONG for a view because of the copy/move semantics of view.
// Use swap instead (the correct one, found by ADL).
namespace std { 
 template <typename V, int R,  typename Opt >
  void swap( triqs::arrays::array_view<V,R,Opt> & a , triqs::arrays::array_view<V,R,Opt> & b)= delete;
}

#ifdef TRIQS_HAS_LAZY_EXPRESSIONS
#include <triqs/clef/adapters/array.hpp>
#endif

#endif

