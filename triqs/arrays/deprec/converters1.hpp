
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

#ifndef TRIQS_ARRAYS_WITH_PYTHON_SUPPORT
#error "You must define the macro TRIQS_ARRAYS_WITH_PYTHON_SUPPORT to use Python interface"
#endif

//#include <boost/python/slice.hpp>
#include "../python/numpy_interface.hpp"
#include <complex>

/* ***********************************************

   Function to register the convert between numpy object and triqs_arrays into boost.python
Usage : (in the main, before launching the interpreter)
#include <triqs/arrays/python/converters.hpp>"
triqs::arrays::register_boost_converters();

 ****************************************************/

// List all the types for which a converter will be generated.
#define TYPES_FOR_WHICH_CONVERTER_IS_NEEDED (long)(std::complex<double>)(double) 

namespace triqs { namespace arrays { 

 namespace BoostPythonConverters { 

/*
 * namespace python = boost::python;
  using Option::option;

  namespace details {
   inline PyObject* to_pyo( python::object const & obj) { 
    python::object Res =  obj;// keep it until the end of the routine
    PyObject * res = Res.ptr();
    Py_INCREF(res); return res;
   }  
  }





  template<typename T,int rank, typename Opt >
   struct array_to_python { static PyObject* convert (array<T,rank,Opt> const & a) { return numpy_interface::from_array(a); } };

  template<typename T,int rank, typename Opt >
   struct array_view_to_python { static PyObject* convert (array_view<T,rank,Opt> const & a) { return numpy_interface::from_array_view(a); } };

  //-----------------------------------------------------------------------------

  template<typename ArrayOrArrayViewType>
   struct Array_from_python_str {

    Array_from_python_str(){ python::converter::registry::push_back(&convertible, &construct, python::type_id< ArrayOrArrayViewType >()); }

    static void* convertible(PyObject* obj_ptr) {
     try  { ArrayOrArrayViewType A(obj_ptr); }
     catch(std::exception s) { std::cerr<<"Array from python conversion impossible : exception is "<<s.what()<<std::endl; return 0;} 
     return obj_ptr;
    }

    static void construct( PyObject* obj_ptr,
      python::converter::rvalue_from_python_stage1_data* data) {
     typedef python::converter::rvalue_from_python_storage<ArrayOrArrayViewType > storage_t;
     storage_t* the_storage = reinterpret_cast<storage_t*>( data );
     void* memory_chunk = the_storage->storage.bytes;
     new (memory_chunk) ArrayOrArrayViewType(obj_ptr);
     data->convertible = memory_chunk;
    }
   };

  //-----------------------------------------------------------------------------

  template<typename T,int rank, typename Opt >
   void register_converter_python_array() {
    python::to_python_converter< array     <T,rank,Opt>, array_to_python     <T,rank,Opt> >();
    python::to_python_converter< array_view<T,rank,Opt>, array_view_to_python<T,rank,Opt> >();
    Array_from_python_str< array<T,rank,Opt> >();
    Array_from_python_str< array_view<T,rank,Opt> >();
   }

  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------

  struct range_to_python  { 
   static PyObject* convert (range const & a) { return details::to_pyo( python::slice(a.first(), a.last(), a.step())); }  
  };

  //-----------------------------------------------------------------------------

  struct range_from_python_str {

   range_from_python_str(){ python::converter::registry::push_back(&convertible, &construct,python::type_id< range >()); }

   static void* convertible(PyObject* obj_ptr) { return (PySlice_Check(obj_ptr)? obj_ptr : 0); }

   static void construct( PyObject* obj_ptr,
     python::converter::rvalue_from_python_stage1_data* data) {
    typedef python::converter::rvalue_from_python_storage<range > storage_t;
    storage_t* the_storage = reinterpret_cast<storage_t*>( data );
    void* memory_chunk = the_storage->storage.bytes;
    PySliceObject* sl = (PySliceObject*)(obj_ptr);
    new (memory_chunk) range(python::extract<long>(sl->start),python::extract<long>(sl->stop),python::extract<long>(sl->step));
    data->convertible = memory_chunk;
   }
  };
  //-----------------------------------------------------------------------------

  inline void register_converter_python_range() {
   python::to_python_converter< range, range_to_python    >();
   range_from_python_str();
  }
*/

#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>

  template<typename T>
   inline void register_boost_converters_for_one_type() {
#define ONE(z,N,unused)\
    triqs::python_tools::register_converter< triqs::arrays::array<T,BOOST_PP_INC(N),triqs::arrays::Option::C> > ();\
    triqs::python_tools::register_converter< triqs::arrays::array<T,BOOST_PP_INC(N),triqs::arrays::Option::Fortran> > ();\
    triqs::python_tools::register_converter< triqs::arrays::array_view<T,BOOST_PP_INC(N),triqs::arrays::Option::C> > ();\
    triqs::python_tools::register_converter< triqs::arrays::array_view<T,BOOST_PP_INC(N),triqs::arrays::Option::Fortran> > ();
    BOOST_PP_REPEAT(ARRAY_NRANK_MAX, ONE, nil);
#undef ONE
   }

  //-----------------------------------------------------------------------------
 }
 inline void register_boost_converters() {
  //static bool initialized = false;
  //if (initialized) return;
  triqs::python_tools::register_converter< triqs::arrays::range> ();
  //BoostPythonConverters::register_converter_python_range();
#define ONETYPE(r, data, elem) BoostPythonConverters::register_boost_converters_for_one_type<elem>();
  BOOST_PP_SEQ_FOR_EACH(ONETYPE, nil , TYPES_FOR_WHICH_CONVERTER_IS_NEEDED); 
#undef ONETYPE
  //initialized = true;
 }

 }

#undef TYPES_FOR_WHICH_CONVERTER_IS_NEEDED

#endif 

