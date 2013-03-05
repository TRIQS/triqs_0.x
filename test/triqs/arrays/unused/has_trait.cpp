
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

#include <iostream>
#include <vector>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/seq.hpp>

using namespace std;

/**
 * Usage : 
 *  - The macro TRIQS_MAKE_TRAIT_HAS_METHOD ( ReturnType, ArgTypeList, Name, NickName )
 *    defines a trait called has_method_NickName such that : 
 *    has_method_NickName<T>::value is true iff 
 *       ReturnType T::Name(ArgTypeList)
 *    is a valid signature 
 *
 *  - The macro TRIQS_MAKE_TRAIT_HAS_CONST_METHOD ( ReturnType, ArgTypeList, Name, NickName )
 *    defines a trait called has_method_NickName such that : 
 *    has_method_NickName<T>::value is true iff 
 *       ReturnType T::Name(ArgTypeList) const
 *    is a valid signature 
 *
 *  - Where : 
 *     ArgTypeList is a boost preprocessor sequence e.g.
 *     TRIQS_MAKE_TRAIT_HAS_CONST_METHOD (void, (int)(double)(string), Truc, Truc)
 *     makes the trait to check if a object A has method  : 
 *     void A::Truc(int,double,string) const
 *
 */
#define TRIQS_MAKE_TRAIT_HAS_METHOD_AUX(r,s,p,XX) BOOST_PP_COMMA_IF(p) XX   
#define TRIQS_MAKE_TRAIT_HAS_METHOD_IMPL( ReturnType, ArgTypeList, Name, NickName, CONST )\
 template<typename T> struct has_method_##NickName {\
  typedef char yes[1]; typedef char no[2];\
  template<typename U, ReturnType (U::*)(BOOST_PP_SEQ_FOR_EACH_I (TRIQS_MAKE_TRAIT_HAS_METHOD_AUX,~,ArgTypeList)) CONST> struct SFINAE {};\
  template<typename U> static yes& Test(SFINAE<U, &U::Name>*);\
  template<typename U> static no& Test(...);\
  static const bool value = sizeof(Test<T>(0)) == sizeof(yes);\
 };
#define TRIQS_MAKE_TRAIT_HAS_METHOD( ReturnType, ArgTypeList, Name, NickName) TRIQS_MAKE_TRAIT_HAS_METHOD_IMPL(ReturnType,ArgTypeList,Name,NickName,)
#define TRIQS_MAKE_TRAIT_HAS_CONST_METHOD( ReturnType, ArgTypeList, Name, NickName) TRIQS_MAKE_TRAIT_HAS_METHOD_IMPL(ReturnType,ArgTypeList,Name,NickName,const)

//-----------------------------------------

TRIQS_MAKE_TRAIT_HAS_CONST_METHOD(void,(size_t)(size_t)(std::vector<double> &),eval, eval);

/*
 * Usage : 
 * a is an object that may have eval or not
 * synthesize(a) returns an object that has eval : 
 *  - if a has eval, just calls it
 *  - else synthetise it.
 *
 *  Use CRTP !! 
*/
namespace details { 
 template<typename T, bool has_it = has_method_eval<T>::value > struct  synthesize;

 template< typename T>
  struct synthesize<T,true> { 
   T obj;
   synthesize(T const & Obj):obj(Obj){}
   void eval(size_t i, size_t jmax, std::vector<double> & res) const { obj.eval(i,jmax,res); }
  };

 template< typename T> 
  struct synthesize<T,false > { 
   T obj;
   synthesize(T const & Obj):obj(Obj){}
   void eval(size_t i, size_t jmax, std::vector<double> & res) const {
    for (unsigned int u=0; u<jmax; ++u) res[u] = obj(i,u);
   }
  };
}

// optional
//namespace result_of { template< typename T> struct synthesize<T> { details::synthesize<T> type;};}

template<typename T> details::synthesize<T> synthesize( T const & ob) { return details::synthesize<T>(ob);}

//---------------------------------------------
//         EXAMPLE 
//

struct A { 
 double operator()(size_t i, size_t j) const { return (i+10*j);}
};

struct B { 
  void assign_to() { cout<<"ee"<<endl;}
  double operator()(size_t i, size_t j) const { return (i+10*j);}
  void eval(size_t i, size_t jmax, std::vector<double> & res) const {
   for (unsigned int u=0; u<jmax; ++u) res[u] =  2*(*this)(i,u);
  }
};

TRIQS_MAKE_TRAIT_HAS_METHOD(void,,assign_to,assign_to);

#include <algorithm>
#include <iterator>

int main() { 

 cout<<" A ?  "<< has_method_assign_to<A>::value<<endl;
 cout<<" B ?  "<< has_method_assign_to<B>::value<<endl;

 cout<<" A ?  "<< has_method_eval<A>::value<<endl;
 cout<<" B ?  "<< has_method_eval<B>::value<<endl;

 vector<double> res(10);

 A a; B b;
 synthesize(a).eval(1,10,res);
 std::copy(res.begin(), res.end(), std::ostream_iterator<double>(std::cout,"\n"));
 std::cout<<std::endl;

 synthesize(b).eval(1,10,res);
 std::copy(res.begin(), res.end(), std::ostream_iterator<double>(std::cout,"\n"));
 std::cout<<std::endl;


}
