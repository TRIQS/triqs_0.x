
/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011 by M. Ferrero, O. Parcollet
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

#include <boost/python.hpp>
#include <triqs/python_tools/converters.hpp>

using namespace boost::python;
using namespace triqs::python_tools;

#include <iostream>
#include <iterator>
#include <vector>
#include <map>
using namespace std;

void f1( object x) { 

 vector<int> V = Py_to_C::convert< vector<int> >::invoke (x);
 std::ostream_iterator<int> out_it (cout,", "); 
 std::copy(V.begin(), V.end(), out_it);
}

void f11( object x) { 

 vector<vector<int > > V = Py_to_C::convert< vector<vector<int> > >::invoke (x);

 for (size_t i = 0 ; i<V.size(); ++i)
  for (size_t j = 0 ; j<V[i].size(); ++j)
   std::cout<<"--"<<i<<" "<<j<<" "<< V[i][j]<<endl;
}

object f2() { 
 vector<double> V(10,2);
 return C_to_Py::convert< vector<double> >::invoke(V);
}

object f3() { 
 vector< vector< double > > V(3, vector<double>(4,2) );
 return C_to_Py::convert< vector<vector<double> > >::invoke(V);
}

void g(object x) {
 typedef boost::unordered_map<std::string, int > type;
 type M = Py_to_C::convert<type>::invoke(x);
 std::map<std::string,int> res;
 res.insert(M.begin(), M.end());
 for (std::map<std::string,int>::iterator it = res.begin(); it !=res.end(); ++it) {
    std::cout<< it-> first <<" --- > "<< it-> second <<endl;
 }
}

object g2() { 
 typedef boost::unordered_map<std::string, int > type;
 type M; 
 M.insert(std::make_pair("a",10));
 M.insert(std::make_pair("b",28));
 return C_to_Py::convert< type >::invoke(M);
}


void p1( object x) { 
 typedef std::pair< vector<double>, vector<int> > type;
 type res =   Py_to_C::convert< type >::invoke (x);
 std::ostream_iterator<int> out_it (cout,", "); 
 std::copy(res.first.begin(), res.first.end(), out_it);
 cout<<endl;
 std::copy(res.second.begin(), res.second.end(), out_it);
 cout<<endl<<endl;

}


object p2() { 
 vector<double> V(5,2.5);
 vector<int> VV(3,6);
 return C_to_Py::convert <std::pair< vector<double>, vector<int> > >::invoke(std::make_pair(V,VV));
};

typedef std::pair< vector<double>, vector<int> > pair_vd_vi;

pair_vd_vi p3() { 
 vector<double> V(5,4.5);
 vector<int> VV(3,-9);
 return std::make_pair(V,VV);
}


typedef boost::unordered_map<std::string, vector<int > > map_st_vi;

map_st_vi g3() { 
 map_st_vi M; 
 M.insert(std::make_pair("a",vector<int>(3,10)));
 M.insert(std::make_pair("b",vector<int>(2,28)));
 return M; 
}

// test with arrays
#include <triqs/arrays/array.hpp>
#include <triqs/arrays/python/converters.hpp>

using namespace triqs::arrays;
typedef array<double,1> ad1;
typedef std::pair < ad1, ad1 > p_ad1_ad1;

p_ad1_ad1 a1() {
 ad1 a(3), b(2); a()=1; b()=3;
  return std::make_pair ( a,b);
}

BOOST_PYTHON_MODULE(_pytriqs_test_converter) {

 triqs::arrays::register_boost_converters();

 triqs::python_tools::register_converter< pair_vd_vi >();
 triqs::python_tools::register_converter< map_st_vi >();
 triqs::python_tools::register_converter< p_ad1_ad1 >();

 def ("f1",f1);
 def ("f11",f11);
 def ("f2",f2);
 def ("f3",f3);

 def ("g",g);
 def ("g2",g2);
 def ("g3",g3);

 def ("p1",p1);
 def ("p2",p2);
 def ("p3",p3);

 def ("a1",a1);

};
