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
#ifndef TRIQS_ARRAY_IMPL_SLICEABLEOBJECT_H   
#define TRIQS_ARRAY_IMPL_SLICEABLEOBJECT_H

#include "./common.hpp"
#include <boost/type_traits/add_const.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/if.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>
#include <boost/preprocessor/repetition/enum_binary_params.hpp>

namespace triqs { namespace arrays { 

 template <typename ValueType, 
	  typename IndexMapType, 
	  typename Opt, 
	  class ViewTag, 
	  template <class V, int R, class Opt2, class Vtag> class ViewFactory,
	  template<typename IndexMapType2, typename SliceArgs> class Slicer,
	  class Derived > 
	   class sliceable_object   
	   {
	    public : 
	     static const unsigned int rank = IndexMapType::domain_type::rank;

	     /** 
	      * All arguments are simply forwarded to Slicer to produce another indexmap IM2.
	      * If the rank of IM2 is 0, the return type is a (const) reference.
	      * Else it is a view (of the type determined by the tag) with IM2 indexmap, on the same data.
	      */
	     template<bool Const, typename ArgsTuple> struct result_of_call {
	      typedef typename Slicer<IndexMapType,ArgsTuple>::return_type IM2;
	      typedef typename boost::mpl::if_c<Const, typename boost::add_const<ValueType>::type, ValueType>::type V2;
	      static const unsigned int R2 = IM2::domain_type::rank;
	      static const bool rank_is_nul = (R2==0);
	      typedef typename boost::mpl::if_c< rank_is_nul, V2 &, 
		      typename ViewFactory<V2,R2, typename Option::Im_Opt_2_Opt<IM2,Opt>::type, ViewTag >::type >::type type;
	     };

	    private:

	     Derived & self() { return static_cast<Derived & >(*this);}
	     Derived const & self() const { return static_cast<Derived const & >(*this);}

	     template<typename ArgsTuple>      
	      typename boost::enable_if_c  < result_of_call<false, ArgsTuple>::rank_is_nul,
		       typename result_of_call<false,ArgsTuple>::type >::type
			eval_(ArgsTuple const & args) { return self()[args]; } 

	     template<typename ArgsTuple>      
	      typename boost::enable_if_c  < result_of_call<true, ArgsTuple>::rank_is_nul, 
		       typename result_of_call<true,ArgsTuple>::type >::type
			eval_(ArgsTuple const & args) const { return self()[args]; }

	     template<typename ArgsTuple>      
	      typename boost::disable_if_c < result_of_call<false, ArgsTuple>::rank_is_nul, 
		       typename result_of_call<false,ArgsTuple>::type >::type
			eval_(ArgsTuple const & args) { 
			 return typename result_of_call<false,ArgsTuple>::type (
			   Slicer<IndexMapType,ArgsTuple>::invoke(self().indexmap(),args), 
			   self().storage()); 
			}

	     template<typename ArgsTuple>      
	      typename boost::disable_if_c < result_of_call<true, ArgsTuple>::rank_is_nul, 
		       typename result_of_call<true,ArgsTuple>::type >::type
			eval_(ArgsTuple const & args) const { 
			 return typename result_of_call<true,ArgsTuple>::type (
			   Slicer<IndexMapType,ArgsTuple>::invoke(self().indexmap(),args), 
			   self().storage()); 
			}

	    public:

	     /// Equivalent to make_view
	     typename result_of_call<true,boost::tuple<> >::type operator()() const { return self(); } 
	     typename result_of_call<false,boost::tuple<> >::type operator()() { return self(); } 

#ifndef TRIQS_HAS_LAZY_EXPRESSIONS

#ifdef TRIQS_ARRAYS_USE_VARIADIC_TEMPLATES 
	     template<typename... Args>      
	      result_of_call<true,TupleFromVariadic<Args>::type >::type 
	      operator()(Args...  args) const { return eval_( TupleFromVariadic<Args>(args) );}

	     template<typename... Args>      
	      result_of_call<false,TupleFromVariadic<Args>::type >::type 
	      operator()(Args...  args) { return eval_( TupleFromVariadic<Args>(args) );}
#else 
#define IMPL(z, NN, unused)                                \
	     template<BOOST_PP_ENUM_PARAMS(NN, typename T)>\
	     typename result_of_call<true,boost::tuple<BOOST_PP_ENUM_PARAMS(NN, T)> >::type\
	     operator()(BOOST_PP_ENUM_BINARY_PARAMS (NN,T, const & x)) const {\
	      return eval_( boost::make_tuple(BOOST_PP_ENUM_PARAMS(NN, x) ));}\
	     \
	     template<BOOST_PP_ENUM_PARAMS(NN, typename T)>\
	     typename result_of_call<false,boost::tuple<BOOST_PP_ENUM_PARAMS(NN, T)> >::type\
	     operator()(BOOST_PP_ENUM_BINARY_PARAMS (NN,T, const & x)) {\
	      return eval_( boost::make_tuple(BOOST_PP_ENUM_PARAMS(NN, x) ));}
	     BOOST_PP_REPEAT_FROM_TO(1,BOOST_PP_INC(ARRAY_NRANK_MAX) , IMPL, nil);
#undef IMPL
#endif

#else 
	     // Another implementation for the case where we have the lazy expressions ON
#define AUX(z,p,unused) BOOST_PP_COMMA_IF(p) boost::ref(arg##p)
#define IMPL(z, NN, unused)                                \
	     template<BOOST_PP_ENUM_PARAMS(NN, typename T)>\
	     typename boost::disable_if< \
	     boost::mpl::or_<  BOOST_PP_ENUM_BINARY_PARAMS(NN, triqs::lazy::is_lazy<T, > BOOST_PP_INTERCEPT) >, \
	     typename result_of_call<true,boost::tuple<BOOST_PP_ENUM_PARAMS(NN, T)> >::type  >::type\
	     operator()(BOOST_PP_ENUM_BINARY_PARAMS (NN,T, const & x)) const {\
	      return eval_( boost::make_tuple(BOOST_PP_ENUM_PARAMS(NN, x) ));}\
	     \
	     template<BOOST_PP_ENUM_PARAMS(NN, typename T)>\
	     typename boost::disable_if< \
	     boost::mpl::or_<  BOOST_PP_ENUM_BINARY_PARAMS(NN, triqs::lazy::is_lazy<T, > BOOST_PP_INTERCEPT) >, \
	     typename result_of_call<false,boost::tuple<BOOST_PP_ENUM_PARAMS(NN, T)> >::type >::type\
	     operator()(BOOST_PP_ENUM_BINARY_PARAMS (NN,T, const & x)) {\
	      return eval_( boost::make_tuple(BOOST_PP_ENUM_PARAMS(NN, x) ));}

	     BOOST_PP_REPEAT_FROM_TO(1,BOOST_PP_INC(ARRAY_NRANK_MAX) , IMPL, nil);
#undef IMPL
#ifndef NEWNEW	     
#define IMPL(z, NN, unused)                                \
	     template< BOOST_PP_ENUM_PARAMS(NN, typename T) > \
	     typename boost::enable_if< \
	     boost::mpl::or_<  BOOST_PP_ENUM_BINARY_PARAMS(NN, triqs::lazy::is_lazy<T, > BOOST_PP_INTERCEPT) >, \
	     typename boost::proto::result_of::make_expr< boost::proto::tag::function, triqs::lazy::FntDomain , \
	     typename ViewFactory<ValueType,rank, Opt , ViewTag >::type,   \
	     BOOST_PP_ENUM_BINARY_PARAMS(NN,T,const & BOOST_PP_INTERCEPT)>::type  >::type  \
	     operator()( BOOST_PP_ENUM_BINARY_PARAMS(NN,T,const & arg) ) const { \
	      return boost::proto::make_expr<boost::proto::tag::function, triqs::lazy::FntDomain>( \
		typename ViewFactory<ValueType,rank, Opt , ViewTag >::type (self() )  , \
		BOOST_PP_REPEAT(NN,AUX,nil)  );\
	     }

	     BOOST_PP_REPEAT_FROM_TO(1,BOOST_PP_INC(ARRAY_NRANK_MAX) , IMPL, nil);

#undef IMPL
#endif
#undef AUX
#endif
	   };// end class

}}//namespace triqs::arrays 
#endif

