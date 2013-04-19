
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
#include <boost/ref.hpp>
#include "./src/tools/print_typeid.hpp"
using namespace std;

namespace pb { 

 namespace details { 
  template<typename T>
   struct wrap {
    T const & obj;
    wrap( T const & obj_): obj(obj_) {}
    void additional_method() {cout<<"add method"<<endl;}
   };
 }

 template<typename T> details::wrap<T> wrap(T const & x) { return details::wrap<T>(x);}

}

/*
 *  The ok version : makes a copy of the object.
 *  If it is temporary, then it is fine.
 *  If we want to use a ref, use boost::ref.
 */
namespace ok { 
 namespace details { 
  template<typename T>
   struct wrap {
    T obj;
    wrap( T const & obj_): obj(obj_) { cout<<"constructing wrap of type : "<<triqs::arrays::Tools::typeid_name(*this)<<endl;}
    void additional_method() {cout<<"add method"<<endl;}
   };
 }

 template<typename T> details::wrap<T> wrap(T const & x) { return details::wrap<T>(x);}
// template<typename T> details::wrap<boost::reference_wrapper<T> > wrap(boost::reference_wrapper<T> x) { 
//  cout<<"ref wrap call"<<endl; return details::wrap<boost::reference_wrapper<T> >(x);}
}

//---------------------------------------------

struct A { 
 A() {cout<<" construct A"<<endl;}
 A(A const & a) {cout<<" copy construct A"<<endl;}
 ~A() {cout<<" descttruc  A"<<endl;}
};

//---------------------------------------------

using namespace ok;
//using namespace pb;

int main() { 

 wrap(A()).additional_method();
 cout<<" *************************"<<endl;
 {
  details::wrap<A>  w(wrap(A()));

  cout<<" ---"<<endl;

  w.additional_method();
 }
 cout<<" *************************"<<endl;
 {
  A a;
  cout<<" --"<<endl;
  wrap(boost::ref(a)).additional_method();
 }
 cout<<" *************************"<<endl;
 {
  A a;
  cout<<" --"<<endl;
  wrap(a).additional_method();
 }
 cout<<" *************************"<<endl;

}

