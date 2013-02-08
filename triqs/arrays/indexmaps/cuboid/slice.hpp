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
#ifndef TRIQS_ARRAYS_INDEXMAP_CUBOID_SLICE_H
#define TRIQS_ARRAYS_INDEXMAP_CUBOID_SLICE_H
#include <triqs/utility/count_type_occurrence.hpp>
namespace triqs { namespace arrays { namespace indexmaps {  

 typedef size_t         l_type;
 typedef std::ptrdiff_t s_type;

 namespace cuboid_details { 

  template <bool BC> inline void _check_BC ( int N, int ind, size_t B, int ind_min =0) { }
  template <> inline void _check_BC<true> (int N, int ind, size_t B, int ind_min) { 
   if (!((ind >= ind_min) && (ind < int(B)))) TRIQS_ARRAYS_KEY_ERROR << " index "<<N<<" is out of domain: \n    " << ind <<" is not within [0,"<< B <<"[\n";
  }

  template<bool BoundCheck> struct slice_calc { 

   l_type const * li; s_type const * si; l_type * lo; s_type * so; int N, P;

   slice_calc(l_type const * li_,s_type const * si_, l_type * lo_, s_type * so_) : li(li_), si(si_), lo (lo_), so(so_),N(0),P(0) {}

   void one_step(std::ptrdiff_t& offset, size_t R){
    _check_BC <BoundCheck> (N, R, li[N]);
    offset += R*si[N];
   }

   //void one_step(size_t const & li, std::ptrdiff_t const & si, size_t & lo, std::ptrdiff_t & so, std::ptrdiff_t& offset, range R){
   void one_step(std::ptrdiff_t& offset, range R){
    _check_BC<BoundCheck> (N, R.first(),li[N]);
    lo[P]  = ((R.last()==-1 ? li[N] : R.last()) - R.first() + R.step()-1 )/R.step(); // python behaviour
    _check_BC<BoundCheck> (N, R.first() + (lo[P]!=0 ? (lo[P]-1) : 0)*R.step() ,li[N], -1);
    so[P]  = si[N] * R.step();
    offset += R.first() * si[N];
   }

   template<int EllipsisLength, typename Arg0, typename...  Args>
    typename std::enable_if<((EllipsisLength>0) && std::is_base_of<ellipsis, Arg0 >::type::value), void>::type 
    invoke(s_type & offset, Arg0 const & arg0, Args const & ... args ) {
     static constexpr bool dP  = (std::is_base_of<range,Arg0>::type::value ? 1 : 0); // Arg0 is range or ellipsis
     one_step(offset, arg0);
     N++; P +=dP;
     invoke<EllipsisLength-1>(offset, arg0, args...);
    }

   template<int EllipsisLength, typename Arg0, typename...  Args>
    typename std::enable_if<!((EllipsisLength>0) && std::is_base_of<ellipsis, Arg0 >::type::value), void>::type 
    invoke(s_type & offset, Arg0 const & arg0, Args const & ... args ) {
     static constexpr bool dP  = (std::is_base_of<range,Arg0>::type::value ? 1 : 0); // Arg0 is range or ellipsis
     one_step(offset, arg0);
     N++; P +=dP;
     invoke<EllipsisLength>(offset, args...);
    }

   template<int EllipsisLength> void invoke(s_type & offset ) {}

   /*template<int EllipsisLength, typename Arg0, typename...  Args>
     static typename std::enable_if<((EllipsisLength>0) && std::is_base_of<ellipsis, Arg0 >::type::value), void>::type 
     invoke(l_type const * li,s_type const * si, l_type * lo, s_type * so, s_type & offset, Arg0 const & arg0, Args const & ... args ) {
     static constexpr bool dP  = (std::is_base_of<range,Arg0>::type::value ? 1 : 0); // Arg0 is range or ellipsis
     one_step(*li,*si,*lo,*so, offset, arg0);
   //l_type const * li2 =li +1;
   slice_calc<N+1,BoundCheck>::invoke<EllipsisLength-1>(li,si,lo,so, offset, arg0, args...);
   //slice_calc<N+1,BoundCheck>::invoke<EllipsisLength-1>(li+1,si+1,lo+dP,so+dP, offset, arg0, args...);
   }

   template<int EllipsisLength, typename Arg0, typename...  Args>
   static typename std::enable_if<!((EllipsisLength>0) && std::is_base_of<ellipsis, Arg0 >::type::value), void>::type 
   invoke(l_type const * li,s_type const * si, l_type * lo, s_type * so, s_type & offset, Arg0 const & arg0, Args const & ... args ) {
   static constexpr bool dP  = (std::is_base_of<range,Arg0>::type::value ? 1 : 0); // Arg0 is range or ellipsis
   one_step(*li,*si,*lo,*so, offset, arg0);
   slice_calc<N+1,BoundCheck>::invoke<EllipsisLength>(li+1,si+1,lo+dP,so+dP, offset, args...);
   }

   template<int EllipsisLength>
   static void invoke(l_type const * li,s_type const * si, l_type * lo, s_type * so, s_type & offset ) {}
   */

  };

  }//namespace cuboid_details

  // special case of no argument : 
  template<int R, bool BC> struct slicer < cuboid::map<R, BC> >  { typedef cuboid::map < R, BC > r_type; }; 

  // general case
  template<int R,  bool BC, typename... Args> struct slicer < cuboid::map<R, BC>,  Args...>  { 

   static const unsigned int len = sizeof...(Args);
   static_assert((count_type_occurrence<ellipsis,Args...>::value < 2), "Only one ellipsis is permitted");
   static_assert((len>=R || (count_type_occurrence<ellipsis,Args...>::value > 0)), "Too few arguments in slice");
   static_assert(len<=R, "Too many arguments in slice");
   typedef cuboid::map < count_type_occurrence<range,Args...>::value , BC > r_type; 
   //typedef cuboid_map < typename index_order::sliced_memory_order<IO,Args...>::type, BC > r_type; 

   static r_type invoke (cuboid::map<R, BC> const & X, Args ... args) { 
    mini_vector<l_type,r_type::rank> newlengths;
    mini_vector<s_type,r_type::rank> newstrides;
    s_type newstart= X.start_shift();
    constexpr int EllipsisLength = R - len;
    cuboid_details::slice_calc<BC>(&X.lengths()[0],&X.strides()[0],&newlengths[0],&newstrides[0]).template invoke<EllipsisLength>(newstart, args...);
    return r_type(newlengths,newstrides,newstart);// use move construction ?
   };
  }; 

 }}}//namespaces triqs::arrays::indexmaps
#endif
