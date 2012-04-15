
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

#ifndef TRIQS_ARRAY_OPTIONS_H
#define TRIQS_ARRAY_OPTIONS_H
#include "../indexmaps/cuboid/cuboid_map.hpp"
#include "../indexmaps/cuboid/cuboid_slice.hpp"

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <boost/mpl/map.hpp>
#include <boost/mpl/has_key.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/at.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/arithmetic/dec.hpp> 

namespace triqs { namespace arrays {

 namespace Tag { struct NoBoundCheck{}; struct BoundCheck{}; }

 /*
  * Option::option<T1, T2, T3 >
  *
  * Parameters can be given in any order 
  *   Order in memory : Tag::C [default], Tag::Fortran, memory_order<0,2,1>, memory_order_p < a_permutation >
  *   Storage : Tag::shared_block [default]
  *   BoundHandler : NoBoundCheck [default], BoundCheck
  */

 namespace Option { 

  namespace mpl= boost::mpl;

  struct _is_cmo {};
#define TEXT(z, nval, text) int n##nval=-1
#define TEXT2(z, nval, text) n##nval
  template < BOOST_PP_ENUM ( BOOST_PP_DEC(ARRAY_NRANK_MAX), TEXT, dd) >  
   struct memory_order : _is_cmo { typedef Permutations::permutation< BOOST_PP_ENUM ( BOOST_PP_DEC(ARRAY_NRANK_MAX), TEXT2, dd) > perm; }; 
#undef TEXT
#undef TEXT2

  template< class P> struct memory_order_p : _is_cmo { typedef P perm; };

  struct _OrderTag {}; struct _StorageTag {}; struct _BoundTag {}; struct _InitTag{};

  typedef mpl::map<
   mpl::pair<void, void> 
   , mpl::pair< Tag::C, _OrderTag >
   , mpl::pair< Tag::Fortran, _OrderTag >
   , mpl::pair< Tag::shared_block, _StorageTag >
   , mpl::pair< Tag::NoBoundCheck, _BoundTag >
   , mpl::pair< Tag::BoundCheck, _BoundTag >
   , mpl::pair< Tag::no_init, _InitTag >
   , mpl::pair< Tag::nan_inf_init, _InitTag >
   , mpl::pair< Tag::default_init, _InitTag >
   > params_list;

  template <class T, class Enable=void> struct recognize;

  template <class T> struct recognize<T, typename boost::disable_if<boost::is_base_of<_is_cmo,T> >::type > { 
   static_assert( (mpl::has_key<params_list, T>::type::value) , " Parameter unknown");
   typedef typename mpl::pair < typename mpl::at<params_list, T>::type, T > type;
  };

  template <class T> struct recognize<T, typename boost::enable_if<boost::is_base_of<_is_cmo,T> >::type > { 
   typedef typename mpl::pair < _OrderTag, typename T::perm >  type;
  };

  template<class s, class key, class def> 
   struct at2 { 
    typedef typename mpl::at<s,key>::type T;
    typedef typename boost::mpl::if_< boost::is_same<T,mpl::void_> , def, T>::type type;
   };

  template <class T1=void, class T2=void, class T3=void, class T4=void> 
   struct options { 

    typedef mpl::map< 
     typename recognize<T1>::type 
     , typename recognize<T2>::type 
     , typename recognize<T3>::type 
     , typename recognize<T4>::type 
     > opt_list;

    static_assert ( (mpl::count<opt_list,_OrderTag>::type::value <= 1), "Too many OrderTag parameters in Option");
    static_assert ( (mpl::count<opt_list,_StorageTag>::type::value <= 1), "Too many StorageTag parameters in Option");
    static_assert ( (mpl::count<opt_list,_BoundTag>::type::value <= 1), "Too many BoundHandler parameters in Option");

    typedef typename at2<opt_list, _OrderTag, Tag::C>::type IndexOrderTag;
    typedef typename at2<opt_list, _StorageTag, Tag::shared_block>::type StorageTag;

#ifdef TRIQS_ARRAYS_ENFORCE_BOUNDCHECK
    typedef typename at2<opt_list, _BoundTag, Tag::BoundCheck >::type BoundHandlerTag;
#else
    typedef typename at2<opt_list, _BoundTag, Tag::NoBoundCheck >::type BoundHandlerTag;
#endif
    static const bool BoundCheck = boost::is_same<BoundHandlerTag, Tag::BoundCheck>::value;

#ifdef TRIQS_ARRAYS_ENFORCE_INIT_NAN_INF
    typedef typename at2<opt_list, _InitTag, Tag::nan_inf_init >::type InitTag;
#else
#ifdef TRIQS_ARRAYS_ENFORCE_INIT_DEFAULT
    typedef typename at2<opt_list, _InitTag, Tag::default_init >::type InitTag;
#else
    typedef typename at2<opt_list, _InitTag, Tag::no_init >::type InitTag;
#endif
#endif
   };


  typedef options<Tag::C> C;
  typedef options<Tag::Fortran> Fortran;
  typedef C Default;

  using indexmaps::cuboid_map; namespace IndexOrder = indexmaps::IndexOrder;

  // compute the indexmap from the dimension and the IndexOrderTag.
  template<int D, typename P, bool BC> struct dim_iotag_2_imap { 
   static_assert( (Tag::check<Permutations::perm_tag,P>::value), "IndexOrderTag must be Tag::C, Tag::Fortran or a Permutation");
   static_assert( (Permutations::length<P>::value==D), "Error : dimension and size of the permutation must be equal");
   typedef cuboid_map<IndexOrder::custom<P>, BC > type;
  };
  template<int D, bool BC> struct dim_iotag_2_imap<D,Tag::C,BC> { typedef cuboid_map<IndexOrder::C<D>, BC >  type;};
  //template<int D> struct dim_iotag_2_imap<D,use_default,BC> { typedef cuboid_map<IndexOrder::C<D> >  type;};
  template<int D, bool BC> struct dim_iotag_2_imap<D,Tag::Fortran,BC> { typedef cuboid_map<IndexOrder::Fortran<D>, BC >  type;};

  // given Rank and Opt, return the IndexMap
  template<int R, class Opt>  struct R_Opt_2_IM : dim_iotag_2_imap< R, typename Opt::IndexOrderTag, Opt::BoundCheck > {};

  //overrule for Fortran<1> -> C<1> // broken
  //template<class Opt>  struct R_Opt_2_IM<1,Opt> : dim_iotag_2_imap< 1, Tag::C, Opt::BoundCheck > {};

  // inverse function : from an IM, an Opt, compute the new Opt such that it would give this indexmap
  template<class IM> struct IM_2_Ordertag;
  template<int R, bool BC> struct IM_2_Ordertag < cuboid_map<IndexOrder::C<R>, BC > > { typedef Tag::C type;};
  template<int R, bool BC> struct IM_2_Ordertag < cuboid_map<IndexOrder::Fortran<R>, BC > > { typedef Tag::Fortran type;};
  template<class P, bool BC> struct IM_2_Ordertag < cuboid_map<IndexOrder::custom<P>, BC > > { typedef memory_order_p<P> type;};

  template<class IM, class Opt1>  struct Im_Opt_2_Opt { 
   typedef options< typename IM_2_Ordertag<IM>::type, typename Opt1::StorageTag, typename Opt1::BoundHandlerTag> type;
  };

 }// Option namespace

 using Option::R_Opt_2_IM; //using Option::Im_Opt_2_Opt;

 template < class V, int R, class Opt, class ViewTag > struct ViewFactory;

}}//namespace triqs::arrays 
#endif

