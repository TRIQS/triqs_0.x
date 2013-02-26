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
#ifndef TRIQS_ARRAYS_BOOST_CONVERTER_H
#define TRIQS_ARRAYS_BOOST_CONVERTER_H
#include "./array_view_to_python.hpp"
#include <boost/python.hpp>
#include <triqs/python_tools/converters.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include "triqs/python_tools/converters.hpp"
#include <boost/python/handle.hpp>
#include <boost/python/slice.hpp>
#include "../array.hpp"
#include "../matrix.hpp"
#include "../vector.hpp"

namespace triqs { namespace python_tools { 

  template<> struct converter<arrays::range> { 
    static bool Py2C_is_possible (bpy::object obj) { return PySlice_Check(obj.ptr());} 
    static arrays::range Py2C(bpy::object obj) { 
     using bpy::extract;
     PySliceObject* sl = (PySliceObject*)(obj.ptr());     
     return arrays::range(extract<long>(sl->start),extract<long>(sl->stop),extract<long>(sl->step));
    }
    static bpy::object C2Py (arrays::range const & a) { return bpy::slice(a.first(), a.last(), a.step());}
   };

  template<typename T, int rank, ull_t Opt, ull_t To > 
   struct converter<arrays::array_view<T,rank,Opt,To>  > {
    static bool Py2C_is_possible (bpy::object obj) {
     try  { arrays::array_view<T,rank,Opt,To>  A(obj.ptr()); }
     catch(runtime_error s) { std::cerr<<s.what()<<std::endl; return false;}
     return true;
    }
    static arrays::array_view<T,rank,Opt,To> Py2C(bpy::object X) { return  arrays::array_view<T,rank,Opt,To> (X.ptr());}
    static bpy::object C2Py (arrays::array_view<T,rank,Opt,To>  const & a) { 
     PyObject * p = arrays::numpy_interface::array_view_to_python(a);
     return bpy::object( bpy::handle<>(p)); // p is a NEW reference ...
    }
   };

  template<typename T, ull_t Opt, ull_t To > 
   struct converter<arrays::matrix_view<T,Opt,To>  > {
    static bool Py2C_is_possible (bpy::object obj) {
     try  { arrays::matrix_view<T,Opt,To>  A(obj.ptr()); }
     catch(runtime_error s) { return false;}
     return true;
    }
    static arrays::matrix_view<T,Opt,To> Py2C(bpy::object X) { return  arrays::matrix_view<T,Opt,To> (X.ptr());}
    static bpy::object C2Py (arrays::matrix_view<T,Opt,To> a){return converter<arrays::array_view<T,2,Opt,To> >::C2Py(a);} 
   };

  template<typename T, ull_t Opt> 
   struct converter<arrays::vector_view<T,Opt>  > {
    static bool Py2C_is_possible (bpy::object obj) {
     try  { arrays::vector_view<T,Opt>  A(obj.ptr()); }
     catch(runtime_error s) { return false;}
     return true;
    }
    static arrays::vector_view<T,Opt> Py2C(bpy::object X) { return  arrays::vector_view<T,Opt> (X.ptr());}
    static bpy::object C2Py (arrays::vector_view<T,Opt> a){return converter<arrays::array_view<T,1,Opt> >::C2Py(a);} 
   };

  template<typename T, int rank, ull_t Opt, ull_t To > 
   struct converter<arrays::array<T,rank,Opt,To>  > {
    typedef void value_not_suitable_for_automatic_conversion_tag; 
    static bool Py2C_is_possible (bpy::object obj) { 
     try  { arrays::array<T,rank,Opt,To>  A(obj.ptr()); }
     catch(runtime_error s) { std::cerr<<s.what()<<std::endl; return false;}
     return true;
    }
    /*static bool Py2C_is_possible (bpy::object obj) { 
      std::cerr<< "--------------------------------------------------------\n"
      << "Error : trying to convert a Python object into an array " << std::endl
      << "  This is not possible. Only array_view is convertible" << std::endl
      << "--------------------------------------------------------\n";
      return false;} 
      */
    static arrays::array<T,rank,Opt,To> Py2C(bpy::object X) { return arrays::array<T,rank,Opt,To> (X.ptr());}
    static bpy::object C2Py (arrays::array<T,rank,Opt,To>  const & a) { return converter<arrays::array_view<T,rank,Opt,To> >::C2Py(a);}
   };

  template<typename T, ull_t Opt, ull_t To > 
   struct converter<arrays::matrix<T,Opt,To>  > {
    typedef void value_not_suitable_for_automatic_conversion_tag; 
    static bool Py2C_is_possible (bpy::object obj) { 
     try  { arrays::matrix<T,Opt,To>  A(obj.ptr()); }
     catch(runtime_error s) { std::cerr<<s.what()<<std::endl; return false;}
     return true;
    }
    /*static bool Py2C_is_possible (bpy::object obj) { 
      std::cerr<< "--------------------------------------------------------\n"
      << "Error : trying to convert a Python object into a matrix " << std::endl
      << "  This is not possible. Only matrix_view is convertible" << std::endl
      << "--------------------------------------------------------\n";
      return false;} 
      */
    static arrays::matrix<T,Opt,To> Py2C(bpy::object X) { return arrays::matrix<T,Opt,To> (X.ptr());}
    static bpy::object C2Py (arrays::matrix<T,Opt,To> const & a ){return converter<arrays::matrix_view<T,Opt,To> >::C2Py(a);} 
   };

  template<typename T, ull_t Opt> 
   struct converter<arrays::vector<T,Opt>  > {
    typedef void value_not_suitable_for_automatic_conversion_tag; 
    static bool Py2C_is_possible (bpy::object obj) { 
     try  { arrays::vector<T,Opt>  A(obj.ptr()); }
     catch(runtime_error s) { std::cerr<<s.what()<<std::endl; return false;}
     return true;
    }
    /*static bool Py2C_is_possible (bpy::object obj) {
     std::cerr<< "--------------------------------------------------------\n"
      << "Error : trying to convert a Python object into a vector " << std::endl
      << "  This is not possible. Only vector_view is convertible" << std::endl
      << "--------------------------------------------------------\n";
      return false;}
      */
    static arrays::vector<T,Opt> Py2C(bpy::object X) { return arrays::vector<T,Opt> (X.ptr());}
    static bpy::object C2Py (arrays::vector<T,Opt> const & a){return converter<arrays::vector_view<T,Opt> >::C2Py(a);} 
   };

}// python_tools

/* ***********************************************
   Mainly for debug : in real code, it is better to list explicitely the converters needed.
   Function to register the convert between numpy object and triqs_arrays into boost.python
Usage : (in the main, before launching the interpreter)
#include <triqs/arrays/python/converters.hpp>"
arrays::register_boost_converters();

 ****************************************************/
// List all the types for which a converter will be generated.
#define TYPES_FOR_WHICH_CONVERTER_IS_NEEDED (long)(std::complex<double>)(double) 

namespace arrays { 

 template<typename T>
  inline void register_boost_converters_for_one_type() {
#define ONE(z,N,unused)\
   python_tools::register_converter< arrays::array_view<T,BOOST_PP_INC(N)> > ();\
   python_tools::register_converter< arrays::array_view<T,BOOST_PP_INC(N)> > ();
   BOOST_PP_REPEAT(ARRAY_NRANK_MAX, ONE, nil);
#undef ONE

#define CONVERT_ALSO_MATRIX_VECTOR 
#ifdef CONVERT_ALSO_MATRIX_VECTOR
   python_tools::register_converter< arrays::matrix_view<T> > ();
   python_tools::register_converter< arrays::vector_view<T> > ();
#endif
  }

 inline void register_boost_converters() {
  python_tools::register_converter< arrays::range> ();
#define ONETYPE(r, data, elem) register_boost_converters_for_one_type<elem>();
  BOOST_PP_SEQ_FOR_EACH(ONETYPE, nil , TYPES_FOR_WHICH_CONVERTER_IS_NEEDED); 
#undef ONETYPE
 }
}}
#undef TYPES_FOR_WHICH_CONVERTER_IS_NEEDED
#endif 

