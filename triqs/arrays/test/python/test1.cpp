
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
#include <Python.h>
#include "./src/array.hpp"
#include <iostream>
#include <boost/python.hpp>
#include "./src/python/converters.hpp"

using boost::python::object;
using namespace std;
using namespace triqs::arrays;
namespace tqa = triqs::arrays;

template<typename T> void f(T const & x) { std::cout<<" x = "<<x<<endl;}

template<> void f(object const & x) { std::cout<<" python object" <<endl;}

class AA { int x; };

template<> void f(AA const & x) { std::cout<<" AA object" <<endl;}


template<class T>
boost::python::object print_array_view( array_view<T,2 > M ) { 
 cout<<" The array is "<<M<<endl;
 return boost::python::object();
}

template<class T>
boost::python::object print_array( triqs::arrays::array<T,2 > M ) { 
 cout<<" The array is "<<M<<endl;
 return boost::python::object();
}

//--------------------------------

//boost::python::object make_array( ) { 
PyObject * make_array( ) { 

 PyObject * r =  NULL; //Py_None;

 {
 triqs::arrays::array<int,2 > M (2,3); 
for (int i =0; i<2; ++i)
  for (int j=0; j<3; ++j) 
  M(i,j) = 10*i +j + 7;

 cout<<" The array is "<<M<<endl;
 array_view<int,2> V(M);

 r = numpy_interface::array_view_to_python (V);

 cerr<< " r ref "<< r->ob_refcnt<<endl;

  //array_view<int,2> A2 = array_view_from_numpy< array_view<int,2> > (r);
  array_view<int,2> A2 (r); 
cerr << A2<< endl;
 cerr<< " r ref "<< r->ob_refcnt<<endl;
}

cerr<< " fini r ref "<< r->ob_refcnt<<endl;

return r;

}


//-----------------------------------

boost::python::object test1( boost::python::object  M ) { 

 try { 
 cerr<<" start of test1 "<< M.ptr()->ob_refcnt<<endl;

 array_view<long,2> A2(M.ptr());
 array_view<long,2> V(A2);
 cout << V<< endl;

 A2(0,0) = 100;

 cout<<" array python "<<A2<<endl;
 }

 catch (std::string s) { cout<<" Exception (string) ! : "<<s <<endl;}
 catch (const char *c) { cout<<" Exception ! : "<<std::string(c)<<endl;}
 cerr<<" end of test1 "<< M.ptr()->ob_refcnt<<endl;
 return M;
}


boost::python::object test2( boost::python::object  M, tqa::range const & R1, tqa::range const & R2) { 
 cerr<<" entring test2 "<<endl;
 array_view<long,2> A (M.ptr());
 array_view<long,2> V(A(R1,R2));
 cerr<<" array python "<<A<<endl<<R1<<" return view  "<<V<<endl;
 return boost::python::object(V);
}

boost::python::object test3(tqa::range const & R) { return boost::python::object(R);}

boost::python::object test4( boost::python::object  M, tqa::range const & R1, int j) { 
 cerr<<" entring test4 "<<endl;
 array_view<long,2> A(M.ptr());
 array_view<long,1> V(A(R1,j));
 cerr<<" array python "<<A<<endl<<" return view  "<<V<<endl;
 return boost::python::object(V);
}

boost::python::object test5( boost::python::object  M, int i, tqa::range const & R1) { 
 array_view<long,2> A(M.ptr());
 array_view<long,1> V(A(i,R1));
 return boost::python::object(V);
}

boost::python::object test6() { 
 triqs::arrays::array<long,2> A(2,3);
 for (int i =0; i<2; ++i)
  for (int j=0; j<3; ++j) 
   A(i,j) = 10*i+ j;

 triqs::arrays::array_view<long,2> V(A(tqa::range(0,2),tqa::range(0,1)));
 triqs::arrays::array_view<long,2> V0(A);
 //cout<<A<<V<<endl;
 return boost::python::make_tuple(triqs::python_tools::make_object(A),V);
}

#include <boost/python.hpp>
#include <boost/python/def.hpp>

BOOST_PYTHON_MODULE(_array_tests)
{
 triqs::python_tools::register_converter< tqa::range >();
 triqs::python_tools::register_converter< tqa::array_view<long,2> >();
 triqs::python_tools::register_converter< tqa::array_view<long,1> >();
 triqs::python_tools::register_converter< tqa::array_view<double,2> >();
 triqs::python_tools::register_converter< tqa::array_view<double,1> >();
// triqs::python_tools::register_converter< tqa::array<long,2> >();
// triqs::python_tools::register_converter< tqa::array<long,1> >();
// triqs::python_tools::register_converter< tqa::array<double,2> >();
// triqs::python_tools::register_converter< tqa::array<double,1> >();
 using namespace boost::python;
 def("test1", test1, "DOC");
 def("test2", test2, "DOC");
 def("test3", test3, "DOC");
 def("test4", test4, "DOC");
 def("test5", test5, "DOC");
 def("test6", test6, "DOC");
 def("print_array_i", print_array_view<long>, "DOC");
 def("print_array_f", print_array_view<double>, "DOC");
 //def("print_array_i", print_array<long>, "DOC");
 //def("print_array_f", print_array<double>, "DOC");
 def("print_array_view_i", print_array_view<long>, "DOC");
 def("print_array_view_f", print_array_view<double>, "DOC");
 def("make_array", make_array , "DOC");

 def("F", f<int>);
 def("F", f<double>);
 def("F", f<AA>);
 //def("F", f<object>);

 class_<AA>("AA", init<>()) ;

}

