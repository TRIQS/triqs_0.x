/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011-2013 by O. Parcollet
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
#ifndef TRIQS_ARRAYS_INDEXMAP_CUBOID_DOMAIN_H
#define TRIQS_ARRAYS_INDEXMAP_CUBOID_DOMAIN_H
#include "../common.hpp"
#include "../range.hpp"
#include "../permutation.hpp"
#include "../../impl/mini_vector.hpp"
#include "../../impl/exceptions.hpp"
#include <iostream>
#include <sstream>

namespace triqs { namespace arrays { namespace indexmaps { namespace cuboid {
 using namespace triqs::arrays::permutations;//::identity;
 

 /// Standard hyper_rectangular domain for arrays
 template<int Rank>
  class domain_t {
   typedef mini_vector<size_t,Rank> n_uple;
   n_uple lengths_;
   friend class boost::serialization::access;
   template<class Archive> void serialize(Archive & ar, const unsigned int version) { ar & boost::serialization::make_nvp("dimensions",lengths_);}
   public :
   //static const unsigned int rank = Rank;
   static constexpr int rank = Rank;
   typedef n_uple index_value_type;
   domain_t ():lengths_(){}
   domain_t (n_uple const & lengths):lengths_(lengths) {}
   domain_t (mini_vector<int,Rank> const & lengths):lengths_(lengths) {}
   domain_t (std::vector<std::size_t> const & l):lengths_() {
    if (!l.size()==rank) TRIQS_RUNTIME_ERROR << "cuboid domain_t construction : vector size incorrect : got "<<l.size() <<" while expected "<< rank;
    lengths_ = n_uple(l);
   }
   domain_t (const domain_t & C):lengths_(C.lengths_){}
   size_t number_of_elements() const { return lengths_.product_of_elements();}
   bool operator==(domain_t const & X) const { return this->lengths_ == X.lengths_;}
   bool operator!=(domain_t const & X) const { return !(this->lengths_ == X.lengths_);}
   n_uple const & lengths() const { return lengths_;}

  /** Generates the value of the indices of a cuboid_domain.  */
   template <ull_t IterationOrder= permutations::identity(Rank) >
    class gal_generator {
      typedef index_value_type indices_type;
      const domain_t * dom;
      indices_type indices_tuple;
      bool atend;
     public:
      gal_generator (const domain_t & P, bool atEnd=false): dom(&P), atend(atEnd) {}
      bool operator==(const gal_generator & IT2) const { assert((IT2.dom == dom)); return ((IT2.atend==atend) );}
      bool operator!=(const gal_generator & IT2) const { return (!operator==(IT2));}
      indices_type const & operator *() const { return indices_tuple;}
      operator bool () const { return !atend;}
      gal_generator & operator++(){ assert(!atend); atend = advance_impl(std::integral_constant<int,0>()); return *this;}
     private:
      template<int r> bool advance_impl(std::integral_constant<int,r>) {
      constexpr int p = permutations::apply(IterationOrder, r);
      if (indices_tuple[p] < dom->lengths()[p]-1) { ++(indices_tuple[p]); return false;}
      indices_tuple[p] = 0;
      return advance_impl(std::integral_constant<int,r+1>());
      }
      bool advance_impl(std::integral_constant<int,rank>) { return true;}
    };

   typedef gal_generator<> generator;

   generator begin() const { return generator(*this,false);}
   generator end() const   { return generator(*this,true);}
   /* End of generator */

   // Check that key in in the domain
   template<class KeyType> void assert_key_in_domain(KeyType const & key) const {
    std::stringstream fs;
    bool res = key_check_impl(std::integral_constant<int,0>(), key,this->lengths_,fs);
    if (!res) TRIQS_ARRAYS_KEY_ERROR << " key out of domain \n" <<fs.str() ;
   }
   template<int r,class KeyType>
    bool key_check_impl (std::integral_constant<int,r>, KeyType const & key, n_uple const & L, std::stringstream & fs ) const {
     using boost::tuples::get;
     bool cond = (  ( size_t(get<r>(key)) < L[r]));
     if (!cond) fs << "key ["<<r<<"] = "<< get<r>(key) <<" is not within [0,"<<L[r]<<"[\n";
     return key_check_impl(std::integral_constant<int,r+1>(), key,L,fs) && cond;
    }
   template<class KeyType> bool key_check_impl (std::integral_constant<int,rank>, KeyType const &, n_uple const &, std::stringstream &) const { return true;}

   // Check that key in in the domain : variadic form. No need for speed optimisation here, it is just for debug
   template<typename ... Args> void assert_key_in_domain_v (Args const & ... args) const { assert_key_in_domain( std::make_tuple(args...));}

   friend std::ostream & operator<<(std::ostream & out, domain_t const & x){return out<<"Cuboid of rank "<<x.rank<<" and dimensions "<<x.lengths();}
  };

 /*/// ------------  tensor product ------------------------
   template<class D1, class D2>
   typename result_of::tensor_product<D1,D2>::type tensor_product( D1 const & d1, D2 const & d2) {
   mini_vector<size_t,D1::rank + D2::rank> res; const int R1(D1::rank), R2(D2::rank);
   for (int u=0;u<R1; ++u) res[u] = d1.lengths()[u];
   for (int u=0;u<R2; ++u) res[D1::rank + u ] = d2.lengths()[u];
   return  domain_t<D1::rank + D2::rank> (res);
   }
   */
}
/// ------------  Pretty Printing : specializing the default behaviour for d=1,2  -------------------------
namespace PrettyPrint_details {
 // TO BE CLEANED
 template<typename A>
  struct print_impl <cuboid::domain_t<1>,A> {
   static void do_it (std::ostream & out,const cuboid::domain_t<1> & d, A const & a ) { out<<"[";
    for (size_t i=0; i< d.lengths()[0]; ++i) out<<(i>0 ? ",": "")<<a(i);
    out<<"]";}
  };

 template<typename A>
  struct print_impl <cuboid::domain_t<2>,A> {
   static void do_it (std::ostream & out,const cuboid::domain_t<2> & d, A const & a ) {
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
}}}
#endif
