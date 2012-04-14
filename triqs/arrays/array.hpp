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
#include "impl/providers.hpp"
#include "impl/option.hpp"

namespace triqs { namespace arrays {

 template <typename ValueType, int Rank, typename Opt= Option::Default > class array_view;
 template <typename ValueType, int Rank, typename Opt= Option::Default > class array;
 
 template <typename ValueType, int Rank, typename Opt>
  class array_view : Tag::array_view, Tag::array_algebra_expression_terminal, 
  public details::indexmap_storage_pair < typename R_Opt_2_IM<Rank,Opt>::type, storages::shared_block<ValueType>, Opt, Tag::array_view >,
  public providers::compound_assign_ops<array_view<ValueType,Rank,Opt> > {
   static_assert( Rank>0, " Rank must be >0");
   public :
    typedef details::indexmap_storage_pair < typename R_Opt_2_IM<Rank,Opt>::type, storages::shared_block<ValueType>, Opt, Tag::array_view > BaseType;
    typedef typename BaseType::indexmap_type indexmap_type;
    typedef array_view<ValueType,Rank,Opt> view_type;
    typedef array<ValueType,Rank,Opt> non_view_type;
    typedef void has_view_type_tag;

    /// Build from an IndexMap and a storage 
    template<typename S> array_view (indexmap_type const & Ind,S const & Mem): BaseType(Ind, Mem) {}

    /// Copy constructor
    array_view(array_view const & X): BaseType(X.indexmap(),X.storage()) {}

    /// Build from anything that has an indexmap and a storage compatible with this class
    template<typename ISP> array_view(const ISP & X): BaseType(X.indexmap(),X.storage()) {}

#ifdef TRIQS_ARRAYS_WITH_PYTHON_SUPPORT
    /**
    * Build from a numpy : only if TRIQS_ARRAYS_WITH_PYTHON_SUPPORT is defined
    */
   explicit array_view (PyObject * X); // implemented in python/numpy_interface : include before use
#endif

    /// Assignment. The size of the array MUST match exactly. 
    template<typename RHS> array_view & operator=(RHS const & X) { assignment(*this,X); return *this; }

    ///
    array_view & operator=(array_view const & X) { assignment(*this,X); return *this; } //without this, the standard = is synthetized...

  };

 template < class V, int R, class Opt > struct ViewFactory< V, R, Opt, Tag::array_view > { typedef array_view<V,R,Opt> type; };

 template <typename ValueType, int Rank, typename Opt>
  class array: Tag::array, Tag::array_algebra_expression_terminal, 
  public  array_view<ValueType,Rank,Opt>::BaseType,
  public providers::compound_assign_ops<array<ValueType,Rank,Opt> >  { 
   typedef typename array_view<ValueType,Rank,Opt>::BaseType BaseType;
   public:
   typedef typename BaseType::value_type value_type;
   typedef typename BaseType::storage_type storage_type;
   typedef typename BaseType::indexmap_type indexmap_type;
   typedef array_view<ValueType,Rank,Opt> view_type;
   typedef array<ValueType,Rank,Opt> non_view_type;
   typedef void has_view_type_tag;

   /// Empty array.
   array():BaseType(indexmap_type(),storage_type()) {}

   /// From a domain
   array( typename indexmap_type::domain_type const & dom):BaseType(indexmap_type(dom)){}

#ifdef TRIQS_DOXYGEN
   /// Construction from the dimensions. NB : the number of parameters must be exactly rank (checked at compile time). 
   array (size_t I_1, .... , size_t I_rank);
#else
#define IMPL(z, NN, unused)                                \
   array (BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(NN), size_t I_)): \
   BaseType(indexmap_type(mini_vector<size_t,BOOST_PP_INC(NN)>(BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(NN), I_)))) {\
    static_assert(BaseType::rank-1==NN,"array : incorrect number of variables in constructor");}
   BOOST_PP_REPEAT(ARRAY_NRANK_MAX , IMPL, nil)
#undef IMPL
#endif

    // Makes a true (deep) copy of the data. 
    array(const array & X): BaseType(X.indexmap(),X.storage().clone()) {}

   /** 
    * Build a new array from X.domain() and fill it with by evaluating X. X can be : 
    *  - another type of array, array_view, matrix,.... (any <IndexMap, Storage> pair)
    *  - a expression : e.g. array<int, IndexOrder::C<2> > A( B+ 2*C);
    */
   template <typename T> 
    array(const T & X, typename boost::enable_if< is_array_assign_lhs<T> >::type *dummy =0):
     BaseType(indexmap_type(X.domain())) { assignment(*this,X); }

#ifdef TRIQS_ARRAYS_WITH_PYTHON_SUPPORT
   /**
    * Build from a numpy : only if TRIQS_ARRAYS_WITH_PYTHON_SUPPORT is defined
    */
   explicit array (PyObject * X); // implemented in python/numpy_interface : include before use
#endif

   /** 
    * Resizes the array. NB : all references to the storage is invalidated.
    * Does not initialize the array by default: to resize and init, do resize(IND).init()
    */
   array & resize (const indexmaps::cuboid_domain<BaseType::rank> & l) { BaseType::resize(l); return *this; }

   /// Assignement resizes the array.  All references to the storage are therefore invalidated.
   array & operator=(const array & X) { BaseType::resize_and_clone_data(X); return *this; }

   /** 
    * Assignement resizes the array (if necessary).
    * All references to the storage are therefore invalidated.
    * NB : to avoid that, do make_view(A) = X instead of A = X
    */
   template<typename RHS> 
    array & operator=(const RHS & X) { 
     static_assert(  is_array_assign_lhs<RHS>::value, "Assignment : RHS not supported");
     BaseType::resize(X.domain());
     assignment(*this,X);
     return *this; 
    }

  };//array class
}}//namespace triqs::arrays

#endif

