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
#ifndef TRIQS_PYTHON_C_CONVERTERS_H 
#define TRIQS_PYTHON_C_CONVERTERS_H
#include <list>
#include <vector>
#include <sstream>
#include <string> 
#include <boost/python/tuple.hpp>
#include <boost/python/list.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/str.hpp>
#include <boost/python/extract.hpp>
#include <boost/unordered_map.hpp>

#include "../utility/typeid_name.hpp"
#include "../utility/exceptions.hpp"

namespace triqs { namespace python_tools { 

 namespace bpy= boost::python;
 inline std::string object_to_string (const bpy::object & O1) { return bpy::extract<std::string>(bpy::str(O1)); }
 inline std::string object_to_string (PyObject * p) { bpy::object obj ( bpy::borrowed (p)); return object_to_string(obj); }
 template<typename T> std::string printer(T const & x) { return triqs::utility::typeid_name(x);}

 // The general case
 template<typename C_Type> struct converter { 

  typedef void _not_suitable_for_automatic_conversion_tag; // avoid infinite loop !

  static bool Py2C_is_possible (bpy::object x) { return bpy::extract<C_Type>(x).check();}

  static C_Type Py2C (bpy::object x) { return  bpy::extract<C_Type> (x);}

  static bpy::object C2Py (C_Type const & x) { 
   bpy::object r;
   try {  r = bpy::object (x); } 
   catch(...) { TRIQS_RUNTIME_ERROR<<"The object "<<printer(x)<<" can not be converted to python "; }
   return r;
  }
 };

 /// More general than boost::object constructor since it does not need registration
 template<typename T> bpy::object make_object(T const & x) { return converter<T>::C2Py(x);}

 //****************** registering the converters in boost python ***********************************8   

 namespace details {
  using namespace boost;

  inline PyObject* to_pyo( python::object const & obj) { 
   python::object Res =  obj;// keep it until the end of the routine
   PyObject * res = Res.ptr();
   Py_INCREF(res); return res;
  }  

  template<typename T, typename Void=void> struct is_suitable_for_automatic_conversion : mpl::true_{};
  template<typename T> struct is_suitable_for_automatic_conversion <T, typename converter<T>::_not_suitable_for_automatic_conversion_tag> : mpl::false_{};
  
  template<typename T, typename Void=void> struct is_suitable_for_automatic_conversion_value : mpl::true_{};
  template<typename T> struct is_suitable_for_automatic_conversion_value <T, typename converter<T>::value_not_suitable_for_automatic_conversion_tag> : mpl::false_{};

  template<class T>
   struct _to_python  { 
    static PyObject* convert (T const & a) { python::object obj =  converter<T>::C2Py(a); return to_pyo(obj); }  
   };

  template<class T> 
   struct _from_python_str {
    _from_python_str(){ python::converter::registry::push_back(&convertible, &construct,python::type_id< T >()); }

    static void* convertible(PyObject* obj_ptr) {
     python::object obj ( python::borrowed (obj_ptr));
     return (converter<T>::Py2C_is_possible(obj) ? obj_ptr : 0);
    } 

    static void construct( PyObject* obj_ptr, python::converter::rvalue_from_python_stage1_data* data) {
     typedef python::converter::rvalue_from_python_storage< T > storage_t;
     storage_t* the_storage = reinterpret_cast<storage_t*>( data );
     void* memory_chunk = the_storage->storage.bytes;
     python::object obj ( python::borrowed (obj_ptr));
     new (memory_chunk) T(converter<T>::Py2C(obj) );
     data->convertible = memory_chunk;
    }
   };

  template<class T>  inline void register_converter() {
   static_assert(is_suitable_for_automatic_conversion<T>::value, "No converter has been found for registering. Did you forget an include ?");
   static_assert(is_suitable_for_automatic_conversion_value<T>::value, "A value (like array, matrix, ...) can not be registered for automatic conversion.");
   static bool initialized = false;
   if (initialized) return;
   bpy::to_python_converter< T, _to_python<T> >(); _from_python_str<T> (); 
   initialized = true;
  }
 } // details

 using details::register_converter;

}}

#endif
