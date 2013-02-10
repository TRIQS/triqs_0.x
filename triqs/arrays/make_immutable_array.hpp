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
#ifndef TRIQS_ARRAYS_MAKE_IMMUTABLE_ARRAY_H
#define TRIQS_ARRAYS_MAKE_IMMUTABLE_ARRAY_H

#ifndef TRIQS_HAS_CLEF_EXPRESSIONS
#error "arrays : make_immutable_xxx can only be used if lazy expressions are included too"
#endif

#include "./array.hpp"
namespace triqs { namespace arrays { 


 // 2 generalization : N dim for array with preproc and matrix/vector.....
 // after rethinking a bit...
 
 template<typename Expr, typename PH1, typename PH2>
  class immutable_array_impl : TRIQS_MODEL_CONCEPT(ImmutableCuboidArray) { 

   public : 
    immutable_array_impl(Expr e_, boost::fusion::pair<PH1,range> p1, boost::fusion::pair<PH2,range> p2): 
     f(clef::make_function(e_, PH1(),PH2())), dom_(make_shape(p1.second.size(), p2.second.size())) {};

    typedef triqs::clef::make_function_impl<Expr, mpl::vector<PH1, PH2> > function_type;
    // pass this result_type in std result of format ...
    typedef typename function_type::template result_type<size_t,size_t>::type value_type; 
    typedef indexmaps::cuboid::domain<2> domain_type;
    domain_type domain() const { return dom_;} 
    template<typename KeyType> value_type operator[] (KeyType const & key) const { return f(key[0],key[1]);}

    friend std::ostream & operator<<(std::ostream & out, immutable_array_impl const & x){return out<<" immutable_array on domain : "<<x.domain();}

   protected:
    function_type f;
    domain_type dom_;

  };

 /** 
  * \brief Makes a lazy immutable (cuboid) array from a simple expression and a set of range...
  * \param expr The lazy expression
  * \param i_=R i_ is a placeholder, R a range. The i_=R produce a (boost::fusion) pair of i_ and R , which is the parameter.
  * \return A lazy object implementing the ImmutableCuboidArray concept with the domain built from the ranges. 
  */
 template<typename Expr, typename PH1, typename PH2>
  immutable_array_impl<Expr,PH1,PH2> make_immutable_array( Expr const & expr, boost::fusion::pair<PH1,range> p1, boost::fusion::pair<PH2,range> p2) {
   return immutable_array_impl<Expr,PH1,PH2> (expr, p1, p2);
  }

}}//namespace triqs::arrays
#endif

