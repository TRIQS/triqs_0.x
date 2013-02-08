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
#ifndef TRIQS_ARRAYS_INDEXMAP_CUBOID_INTERNAL_H
#define TRIQS_ARRAYS_INDEXMAP_CUBOID_INTERNAL_H

namespace triqs { namespace arrays { namespace indexmaps {  

 namespace cuboid_details { 

  template <bool BC> inline void _check_BC ( int N, int ind, size_t B) { }
  template <bool BC> inline void _check_BC2 ( int N, int ind, size_t B) { }

  template <> inline void _check_BC<true> (int N, int ind, size_t B) { 
   bool cond = (ind >= 0) && (ind < int(B));
   if (!cond) TRIQS_ARRAYS_KEY_ERROR << " index "<<N<<" is out of domain: \n    " << ind <<" is not within [0,"<< B <<"[\n";
  }
  template <> inline void _check_BC2<true> (int N, int ind, size_t B) { 
   bool cond = (ind >= -1) && (ind < int(B));
   if (!cond) TRIQS_ARRAYS_KEY_ERROR << " index "<<N<<" is out of domain: \n    " << ind <<" is not within [0,"<< B <<"[\n";
  }

  template<int Rank_in, int Rank_out, int N, int Na, int P, int c, bool BoundCheck> struct slice_calc { 

   typedef mini_vector<size_t,Rank_in> const &  Li_type;
   typedef mini_vector<size_t,Rank_out>  &      Lo_type;
   typedef mini_vector<std::ptrdiff_t,Rank_in > const & Si_type; 
   typedef mini_vector<std::ptrdiff_t,Rank_out> &       So_type; 

   template< typename ArgsTuple>
    static void invoke(Li_type li, Si_type si, Lo_type lo, So_type so, std::ptrdiff_t & offset, ArgsTuple const & args ) {
     one_step(li, si,lo,so, offset,boost::tuples::get<Na>(args));
     static const unsigned int lenA = boost::tuples::length<ArgsTuple>::value;
     static const bool is_range    = boost::is_base_of<typename boost::tuples::element<Na,ArgsTuple>::type, range >::type::value ;
     static const bool is_ellipsis = boost::is_base_of<typename boost::tuples::element<Na,ArgsTuple>::type, ellipsis >::type::value ;
     static const int dP = ( (is_range||is_ellipsis) ? 1 : 0) ;
     // if it is an ellipsis, we do not consume the next argument, unless # of remaining args <= # args to fill, for the right part ...
     static const int dNa = ( ( is_ellipsis && (Rank_in - N  > lenA - Na )) ? 0 : 1) ;
#ifdef TRIQS_ARRAYS_DEBUG_SLICE
     std::cerr<<"N = "<< N << "  P = "<<P<< " c = "<<c<<"dP= "<< dP <<std::endl;
#endif
     slice_calc<Rank_in,Rank_out,N+1,Na + dNa, P+dP, c-1, BoundCheck>::invoke(li,si,lo,so, offset, args);
    }

   static void one_step(Li_type li, Si_type si, Lo_type lo, So_type so,std::ptrdiff_t& offset, size_t R){
    _check_BC <BoundCheck> (N, R, li[N]);
    offset += R*si[N];
   }

   static void one_step(Li_type li, Si_type si, Lo_type lo, So_type so, std::ptrdiff_t& offset, range R){
    _check_BC <BoundCheck> (N, R.first(),li[N]);
    //_check_BC2 <BoundCheck> (N, (R.last()==-1 ? li[N] : R.last()) -1 ,li[N]);
    lo[P]  = ((R.last()==-1 ? li[N] : R.last()) - R.first() + R.step()-1 )/R.step(); // python behaviour
    _check_BC2 <BoundCheck> (N, R.first() + (lo[P]!=0 ? (lo[P]-1) : 0)*R.step() ,li[N]);
    so[P]  = si[N]  * R.step();
    offset += R.first() * si[N]  ;
   }

  };

  // stop the recursion
  template<int Ri, int Ro, int N, int Na, int P, bool BC> struct slice_calc <Ri,Ro,N,Na,P,0,BC> : slice_calc<Ri,Ro,N,Na,P,1,BC> { 
   template<class T1,class T2,class T3,class T4,class T5,class T6> static void invoke(T1,T2,T3,T4,T5,T6 ) {}
  };


 }//namespace cuboid_details

 template<typename IO, typename ArgsTuple, bool BC> 
  struct slicer < cuboid_map<IO, BC>,  ArgsTuple>  { 

   static const unsigned int len = boost::tuples::length<ArgsTuple>::value;
   static_assert((TupleTools::CountHowMany<ellipsis,ArgsTuple>::value < 2), "Only one ellipsis is permitted");
   static_assert((len>=IO::rank || (TupleTools::CountHowMany<ellipsis,ArgsTuple>::value > 0)), "Too few arguments in slice");
   static_assert(len<=IO::rank, "Too many arguments in slice");

   typedef cuboid_map < typename IndexOrder::SlicedIndexOrder<IO, ArgsTuple>::type, BC > return_type; 

   static return_type invoke (cuboid_map<IO, BC> const & X, ArgsTuple args) { 
    mini_vector<size_t        ,return_type::rank> newlengths;
    mini_vector<std::ptrdiff_t,return_type::rank> newstrides;
    std::ptrdiff_t newstart= X.start_shift();
    cuboid_details::slice_calc<IO::rank,return_type::rank,0,0,0,IO::rank, BC>::invoke(X.lengths(),X.strides(),newlengths,newstrides, newstart, args);
#ifdef TRIQS_ARRAYS_DEBUG_SLICE
    std::cerr<<"-----------------------------------------------"<<std::endl;
    std::cerr<<"Slicing "<< X.lengths().to_string()<<X.strides().to_string()<<newlengths.to_string()<<newstrides.to_string()<< newstart << args<<std::endl;
#endif
    return return_type(newlengths,newstrides,newstart);
   };

  }; 

 // special case of no argument : 
 template<typename IO, bool BC> struct slicer < cuboid_map<IO, BC>, boost::tuple<> >  { typedef cuboid_map < IO, BC > return_type; }; 


}}}//namespaces triqs::arrays::indexmaps
#endif
