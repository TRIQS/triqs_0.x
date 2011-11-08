
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

 inline std::string object_to_string (const boost::python::object & O1) {
  std::stringstream fs;
  fs<<std::string(boost::python::extract<std::string>(boost::python::str(O1)));
  return fs.str();
 }

 namespace Py_to_C { 

  template<typename T> struct convert{ 
   static bool possible (boost::python::object x) { boost::python::extract<T> ext (x); return ext.check(); }
   static T invoke(boost::python::object x) { 
    boost::python::extract<T> ext (x);
    if (!ext.check())  {
     TRIQS_RUNTIME_ERROR<<"Object  : "<<object_to_string(x)
     <<"\n of python type "<< object_to_string(x.attr("__class__").attr("__name__"))
     <<"\n can not be converted to the C++ type "<<triqs::utility::typeid_name(T()); 
    }
     return ext();
   }
  };

  // VECTOR
  template<typename T> 
   struct convert<std::vector<T> > { 
    static bool possible (boost::python::object l) { 
     try { 
      boost::python::list L(l); // check it is a list...
      const size_t N = boost::python::len(L);
      for (size_t i=0; i<N; ++i) if (!convert<T>::possible(l[i])) throw false;
      return true;
     }
     catch(...) {return false;}
    }

    static std::vector<T> invoke(boost::python::object l) { 
     boost::python::list L(l);
     const size_t N = boost::python::len(L);
     std::vector<T> res;
     for (size_t i=0; i<N; ++i) res.push_back( convert<T>::invoke(l[i]));
     return res;
    }
   };

  // PAIR
  template<typename T1, typename T2> 
   struct convert<std::pair<T1,T2> > { 

    static bool possible (boost::python::object l) { 
     try { 
      //boost::python::list L(l); // check it is a list...
      if (boost::python::len(l)!=2) throw false;
      if (! (convert<T1>::possible(l[0])) && (convert<T2>::possible(l[1])) ) throw false;
      return true;
     }
     catch(...) {return false;}
    }

    static std::pair<T1,T2> invoke(boost::python::object l) { 
     boost::python::list L(l);
     const size_t N = boost::python::len(L);
     if (N!=2) TRIQS_RUNTIME_ERROR<<"Object"<<object_to_string(l)<<" can not be converted to a std::pair : lengths mismatch "; 
     return std::make_pair( convert<T1>::invoke(l[0]), convert<T2>::invoke(l[1]) );
    }
   };

  // DICT  <-> unordered_map
  template<typename Key, typename Mapped> 
   struct convert < boost::unordered_map<Key,Mapped> >  {

    static bool possible (boost::python::object d) { 
     try { 
      boost::python::dict dic(d); boost::python::list keys = dic.keys(), vals = dic.values();
      const size_t N = boost::python::len(dic);
      for (size_t i=0; i<N; ++i) if ( (!convert<Key>::possible(keys[i])) || (!convert<Mapped>::possible(vals[i]))) throw false;
      return true;
     }
     catch(...) {return false;}
    }

    static boost::unordered_map<Key,Mapped> invoke (boost::python::object d) { 
     boost::unordered_map<Key,Mapped> res;
     boost::python::dict dic(d);
     boost::python::list keys = dic.keys(), vals = dic.values();
     const size_t N = boost::python::len(dic);
     for (size_t i=0; i<N; ++i) {
      bool ok =  res.insert(std::make_pair(convert<Key>::invoke(keys[i]), convert<Mapped>::invoke(vals[i]))).second;
      if (!ok) TRIQS_RUNTIME_ERROR<<"This error should never happen !";
     }
     return res;
    }
   };
 }

 namespace C_to_Py { 

  // a nice printing of the object for debug : type
  template<typename T> std::string printer(T const & x) { return triqs::utility::typeid_name(x);}

  template<typename T> struct convert { 
   static boost::python::object invoke (T const & x) { 
    boost::python::object r;
    try {  r = boost::python::object (x); } 
    catch(...) { TRIQS_RUNTIME_ERROR<<"The object "<<printer(x)<<" can not be converted to python "; }
    return r;
   }
  };

  template<typename T1, typename T2> 
   struct convert<std::pair<T1,T2> > { 
    static boost::python::object invoke (std::pair<T1,T2> const & p) { 
     return boost::python::make_tuple( convert<T1>::invoke (p.first), convert<T2>::invoke (p.second) );
    }
   };

  template<typename T> 
   struct convert< std::vector<T>  >  { 
    static  boost::python::object invoke (std::vector<T> const &v) {
     const size_t N = v.size();
     boost::python::list L;
     for (size_t i=0; i<N; ++i) L.append( convert<T>::invoke (v[i])); 
     return L;
    }
   };


  template<typename Key, typename Mapped> 
   struct convert< boost::unordered_map<Key,Mapped> > { 
    static boost::python::object invoke( boost::unordered_map<Key,Mapped> const & m) {
     boost::python::dict Res;
     for (typename boost::unordered_map<Key,Mapped>::const_iterator it = m.begin(); it !=m.end(); ++it) 
      Res[convert<Key>::invoke (it->first)] = convert<Mapped>::invoke (it->second);
     return Res;
    }
   };

 }


 //****************** registering the converters in boost python ***********************************8   

 namespace details {
  using namespace boost;

  inline PyObject* to_pyo( python::object const & obj) { 
   python::object Res =  obj;// keep it until the end of the routine
   PyObject * res = Res.ptr();
   Py_INCREF(res); return res;
  }  

  template<class T>
   struct _to_python  { 
    static PyObject* convert (T const & a) { python::object obj =  C_to_Py::convert<T>::invoke(a); return to_pyo(obj); }  
   };

  template<class T> 
   struct _from_python_str {
    _from_python_str(){ python::converter::registry::push_back(&convertible, &construct,python::type_id< T >()); }

    static void* convertible(PyObject* obj_ptr) {
     python::object obj ( python::borrowed (obj_ptr));
     return (Py_to_C::convert<T>::possible(obj) ? obj_ptr : 0);
    } 

    static void construct( PyObject* obj_ptr, python::converter::rvalue_from_python_stage1_data* data) {
     typedef python::converter::rvalue_from_python_storage< T > storage_t;
     storage_t* the_storage = reinterpret_cast<storage_t*>( data );
     void* memory_chunk = the_storage->storage.bytes;
     python::object obj ( python::borrowed (obj_ptr));
     new (memory_chunk) T(Py_to_C::convert<T>::invoke(obj) );
     data->convertible = memory_chunk;
    }
   };

  template<class T>  inline void register_converter() {
   static bool initialized = false;
   if (initialized) return;
   //T x; boost::python::extract<T> test(x); 
   //if (test()) // if you try to use it for T= double e.g. ....
   //TRIQS_RUNTIME_ERROR<<"triqs::python_tools : the converter for the type "<<triqs::utility::typeid_name(T())<< " is already active !!"; 
   boost::python::to_python_converter< T, _to_python<T> >(); _from_python_str<T> (); 
   initialized = true;
  }
 } // details

 using details::register_converter;

}}

#endif
