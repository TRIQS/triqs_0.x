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

#ifndef TRIQS_ARRAYS_INDEXMAP_CUBOID_H
#define TRIQS_ARRAYS_INDEXMAP_CUBOID_H
#include "../common.hpp"
#include "../range.hpp"
#include "../../impl/mini_vector.hpp"
#include "../../impl/exceptions.hpp"
#include "./cuboid_index_generator.hpp"
#include <iostream>
#include <sstream>

namespace triqs { namespace arrays { 

 namespace Tag {struct cuboid_domain {};} 
 namespace indexmaps { 
  /// Standard hyper_rectangular domain for arrays
  template<int Rank>
   class cuboid_domain : Tag::cuboid_domain  { 
    typedef mini_vector<size_t,Rank> n_uple;
    public : 
    static const unsigned int rank = Rank;
    typedef n_uple index_value_type;
    cuboid_domain ():lengths_(){}
    cuboid_domain (n_uple const & lengths):lengths_(lengths) { }
    cuboid_domain (mini_vector<int,Rank> const & lengths):lengths_(lengths) { }
    cuboid_domain (std::vector<std::size_t> const & l):lengths_() {
     if (!l.size()==rank) throw "cuboid_domain construction : vector size incorrect";
     lengths_ = n_uple(l); 
    }
    cuboid_domain (const cuboid_domain & C):lengths_(C.lengths_){}
    size_t number_of_elements() const { return lengths_.product_of_elements(); }
    bool operator==(cuboid_domain const & X) const { return this->lengths_ == X.lengths_;}
    bool operator!=(cuboid_domain const & X) const { return !(this->lengths_ == X.lengths_);}
    n_uple const & lengths() const { return lengths_;}
    typedef cuboid_index_generator <cuboid_domain, Permutations::identity<rank> > generator; 
    generator begin() const { return generator(*this,false);}
    generator end() const { return generator(*this,true);}

    template<class KeyType> void assert_key_in_domain(KeyType const & key) const;

    friend std::ostream & operator<<(std::ostream & out, cuboid_domain const & x){return out<<"Cuboid of rank "<<x.rank<<" and dimensions "<<x.lengths(); }

    protected:
    n_uple lengths_;
    friend class boost::serialization::access;
    template<class Archive>
     void serialize(Archive & ar, const unsigned int version) { ar & boost::serialization::make_nvp("dimensions",lengths_); }
   };

  namespace result_of {
   template<class D1, class D2> struct tensor_product { typedef cuboid_domain<D1::rank + D2::rank> type; };
  } 

  /// ------------  tensor product ------------------------
  template<class D1, class D2> 
   typename result_of::tensor_product<D1,D2>::type tensor_product( D1 const & d1, D2 const & d2) { 
    mini_vector<size_t,D1::rank + D2::rank> res; const int R1(D1::rank), R2(D2::rank); 
    for (int u=0;u<R1; ++u) res[u] = d1.lengths()[u];
    for (int u=0;u<R2; ++u) res[D1::rank + u ] = d2.lengths()[u];
    return  cuboid_domain<D1::rank + D2::rank> (res);
   }

  /// ------------  Pretty Printing : specializing the default behaviour for d=1,2  -------------------------
  namespace PrettyPrint_details { 
   template<typename A>
    struct print_impl <cuboid_domain<1>,A> {
     static void do_it (std::ostream & out,const cuboid_domain<1> & d, A const & a ) { out<<"[";
      for (size_t i=0; i< d.lengths()[0]; ++i) out<<(i>0 ? ",": "")<<a(i);
      out<<"]"; }
    };

   template<typename A>
    struct print_impl <cuboid_domain<2>,A> {
     static void do_it (std::ostream & out,const cuboid_domain<2> & d, A const & a ) { 
      out<<"\n[";
      for (size_t i=0; i< d.lengths()[0]; ++i) {
       out<<(i==0 ? "[" : " [");
       for (size_t j=0; j< d.lengths()[1]; ++j) out<<(j>0 ? ",": "")<<a(i,j);
       out<<"]"<<(i==d.lengths()[0]-1 ? "" :  "\n");
      }
      out<<"]";
     }
    };
  }

  namespace cuboid_details {

   template<typename KeyType, typename Ltype, int v> struct key_check_impl {
    static bool invoke ( KeyType const & key, Ltype const & L, std::stringstream & fs ) {
     using boost::tuples::get;
     //bool cond = ( (get<v>(key) >= 0) && ( get<v>(key) < L[v]));
     bool cond = (  ( size_t(get<v>(key)) < L[v]));
     if (!cond) fs << "key ["<<v<<"] = "<< get<v>(key) <<" is not within [0,"<<L[v]<<"[\n";
     return key_check_impl<KeyType,Ltype,v-1>::invoke(key,L,fs) && cond;
    }
   };

   template<typename KeyType, typename Ltype> struct key_check_impl<KeyType,Ltype,-1> {
    static bool invoke ( KeyType const & key, Ltype const & L, std::stringstream & fs ) { return true;}
   };

  }

  template<int Rank>
   template<class KeyType>
   void cuboid_domain<Rank>::assert_key_in_domain(KeyType const & key) const { 
    std::stringstream fs;
    bool res = cuboid_details::key_check_impl<KeyType,n_uple, rank -1>::invoke(key,this->lengths_,fs);
    if (!res) TRIQS_ARRAYS_KEY_ERROR << " key out of domain \n" <<fs.str() ;
   }

 }}}

#endif
