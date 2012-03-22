/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012 by M. Ferrero, O. Parcollet
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

#ifndef TRIQS_TAYLOR_EXPANSION_H
#define TRIQS_TAYLOR_EXPANSION_H 

#include <triqs/utility/algebra.hpp>

namespace triqs { namespace gf { 

 template<typename A>
  class laurent_expansion { 

   const int OrderMinMIN, OrderMaxMAX; // maxi dimension of the expansion. Fixed at construction.

   // add python interface 

   A const & operator[](int i) const { return M[i-OrderMinMIN];}

   protected:
   std::vector<A> M;

  };

 /*
  * ImmutableLaurentExpansion Concept
  *
  * Concept of the formal series \sum_n a_n X^n
  *
  * element_type // type of the a_n
  * element_construct_type // the type to call to construct the element_type (e.g. array vs array_view)
  * ---> ou alors generalize the view_type of the arrays to the triqs:: for everything ! array, GF, etc....
  *
  *  inverse(element_type) // must be accessible by ADL and throw an exception if inversion is not possible.
  *                        
  *  [element_type | element_type const &] operator[int i] const // access to element
  *
  *  const int order_min_min, order_max_max   //
  *  int order_min() const
  *  int order_max() const 
  */

 namespace tag { struct is_laurent_expansion{}; }
  
 template<typename T> struct is_laurent_expansion : boost::is_base_of<tag::is_laurent_expansion,T> {};
 
 // ou avoir un general -> inverse_impl to be specialized ??? better....
 template<typename T> typename boost::enable_if< is_laurent_expansion<T>, T >::type inverse (T const & x);
 
 namespace expressions { 
  
  template<typename T> struct wrap_scalar {
   T s; 
   wrap_neutral( T const & x) : s(x) {}
   typedef T element_type;
   element_type const & operator[](int i) const { return s;}
  };

  template<typename T> struct wrap_negate {
   T const & s; 
   wrap_negate(T const & x) : s(x) {}
   typedef T element_type;
   element_type const & operator[](int i) const { return -s[i];}
  };

  template<typename ProtoTag, typename L, typename R> struct wrap_arith_node : tag::is_laurent_expansion  { 
   L const & l; R const & r;
   wrap_arith_node (L const & l_, R const & r_):l(l_),r(r_) {} 
   typedef triqs::utility::expressions::_ops_<ProtoTag, typename L::element_type, typename R::element_type > ops_type;
   typedef typename ops_type::result_type element_type;
   element_type const & operator[](int i) const { return ops_type::invoke( l[i], r[i]);}
  }; 

  template<typename L, typename R> struct wrap_arith_node<boost::proto::tag::multiplies,L,R> : tag::is_laurent_expansion  { 
   L const & l; R const & r;
   wrap_arith_node (L const & l_, R const & r_):l(l_),r(r_) {} 
   typedef triqs::utility::expressions::_ops_<boost::proto::tag::multiplies, typename L::element_type, typename R::element_type > ops_type;
   typedef typename ops_type::result_type element_type;
   element_type const & operator[](int n) const { 

    const int omin(l.OrderMin()), omax(l.OrderMax());
    const int new_ordermin(omin+r.OrderMin());
    const int new_ordermax(std::min(omax+r.OrderMax(),OrderMaxMAX));

    //if (N1!=t.N1 || N2!=t.N2) TRIQS_RUNTIME_ERROR<<"Multiplication is valid only for similar tail shapes !";
    //if (new_ordermin < OrderMinMIN) TRIQS_RUNTIME_ERROR<<"The multiplication makes the new tail have a too small OrderMin";

    const int kmin(std::max(0,n- r.OrderMax()-omin));
    const int kmax(std::min(omax-omin,n-r.OrderMin()-omin));

    result_type res = typename triqs::non_view_type<typename L::element_type>::type ( l[omin+kmin] * r [n - omin - kmin]);
    for (int k = kmin+1; k <= kmax; ++k) { res += l[omin + k] * r [n - omin -k]; }
    return res;
   }

  };  

  template<typename T> 
   typename boost::enable_if< is_laurent_expansion<T>, std::ostream &  >::type 
   formal_print (std::ostream & out, T const & x) {
    for(size_t u=order_min(); u<=x.order_max(); ++u) out<<x[u] <<" X^"<<u;
   }

  template <typename Expr> struct LaurentExpansionExpr;

  typedef triqs::utility::expressions::grammar_generator<is_laurent_expansion, is_scalar, wrap_scalar, wrap_arith_node, wrap_negate> grammar_generator;
  typedef typename grammar_generator::Grammar LaurentExpansionGrammar;

  struct LaurentExpansionDomain : proto::domain<proto::generator<LaurentExpansionExpr>, LaurentExpansionGrammar> {
   //  template< typename T > struct as_child : proto_base_domain::as_expr< T > {};
  };

  template<typename Expr> struct LaurentExpansionExpr : proto::extends<Expr, LaurentExpansionExpr<Expr>, LaurentExpansionDomain> { 
   typedef proto::extends<Expr, LaurentExpansionExpr<Expr>, LaurentExpansionDomain> base_type;

   typedef typename boost::result_of<LaurentExpansionGrammar(Expr) >::type _T;
   typedef typename _T::element_type element_type;

   LaurentExpansionExpr( Expr const & expr = Expr() ) : base_type( expr ) {}
   value_type operator[] (int n) const { return LaurentExpansionGrammar()(*this)[n]; }
   friend std::ostream &operator <<(std::ostream &sout, LaurentExpansionExpr<Expr> const &expr) { return proto::eval(expr, triqs::utility::expressions::AlgebraPrintCtx (sout)); }
  };
 }

}}

BOOST_PROTO_DEFINE_OPERATORS(triqs::gf::expressions::is_laurent_expansion, triqs::gf::expressions::LaurentExpansionDomain);

#endif
