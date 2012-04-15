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
#ifndef TRIQS_ARRAYS_NUMPY_INTERACE_H
#define TRIQS_ARRAYS_NUMPY_INTERACE_H

#ifdef TRIQS_ARRAYS_WITH_PYTHON_SUPPORT
#include "../array.hpp"
#include "../matrix.hpp"
#include "../vector.hpp"
#include "../storages/shared_block.hpp"
#include "../../utility/exceptions.hpp"
#include "../../utility/typeid_name.hpp"
#include "numpy/arrayobject.h"
#include <boost/python.hpp>

namespace triqs { namespace arrays { namespace numpy_interface  {
 using namespace std;

 inline std::string object_to_string (const boost::python::object & O1) {
  std::stringstream fs;
  fs<<std::string(boost::python::extract<std::string>(boost::python::str(O1)));
  return fs.str();
 }
 inline std::string object_to_string (PyObject * p) { boost::python::object obj ( boost::python::borrowed (p)); return object_to_string(obj); }

 template <class T> struct numpy_to_C_type;
#define CONVERT(C,P) template <> struct numpy_to_C_type<C> { enum {arraytype = P}; }
 CONVERT(bool, NPY_BOOL);
 CONVERT(char, NPY_CHAR);
 CONVERT(signed char, NPY_BYTE);
 CONVERT(unsigned char, NPY_UBYTE);
 CONVERT(short, NPY_SHORT);
 CONVERT(unsigned short, NPY_USHORT);
 CONVERT(int, NPY_INT);
 CONVERT(unsigned int, NPY_UINT);
 CONVERT(long, NPY_LONG);
 CONVERT(unsigned long, NPY_ULONG);
 CONVERT(long long, NPY_LONGLONG);
 CONVERT(unsigned long long, NPY_ULONGLONG);
 CONVERT(float, NPY_FLOAT);
 CONVERT(double, NPY_DOUBLE);
 CONVERT(long double, NPY_LONGDOUBLE);
 CONVERT(std::complex<float>, NPY_CFLOAT);
 CONVERT(std::complex<double>, NPY_CDOUBLE);
 CONVERT(std::complex<long double>, NPY_CLONGDOUBLE);
#undef CONVERT

 struct copy_exception : public triqs::runtime_error {};

 template<typename IndexMapType, typename ValueType >
  std::pair <IndexMapType, storages::shared_block<ValueType> >  extract_len_stri_block( PyObject * X, bool allow_copy) {

   if (X==NULL) TRIQS_RUNTIME_ERROR<<"numpy interface : the python object is NULL !";
   int _r = _import_array();assert(_r==0);

   // make sure IndexMap is cuboid ?s
   static const char Order = IndexMapType::index_order_type::C_or_F;
   static_assert( ((Order=='F')||(Order=='C')), "Ordering must be C or Fortran");

   if ((!PyArray_Check(X)) && (Order=='D'))
    TRIQS_RUNTIME_ERROR<<"numpy interface : the python object is not a numpy and you ask me to deduce the ordering in memory !";

   // From X, we ask the numpy library to make a numpy, and of the correct type.
   // This handles automatically the cases where : 
   //   - we have list, or list of list/tuple
   //   - the numpy type is not the one we want.
   //   - adjust the dimension if needed
   // If X is an array : 
   //   - if Order is same, don't change it
   //   - else impose it (may provoque a copy).
   // if X is not array : 
   //   - Order = FortranOrder or SameOrder - > Fortran order otherwise C
   const int elementsType (numpy_to_C_type<typename boost::remove_const<ValueType>::type>::arraytype);
   int rank = IndexMapType::rank;
   PyObject * numpy_obj;

   if (!allow_copy) {
    // in case of a view, we decide ourselves if we can do it...
    // a previous uses PyArray_FromAny, but behaviour changes between different version of numpy, so it is not portable...
    if (!PyArray_Check(X)) 
     throw copy_exception ()
      <<"   A deep copy of the object would be necessary while views are supposed to guarantee to present a *view* of the python data.\n"
      <<"   Indeed the object was not even an array !\n";  

    if ( elementsType != PyArray_TYPE((PyArrayObject*)X))
     throw copy_exception ()
      <<"   A deep copy of the object would be necessary while views are supposed to guarantee to present a *view* of the python data.\n"
      <<"   The deep copy is caused by a type mismatch of the elements. \n";

    PyArrayObject *arr = (PyArrayObject *)X;    
    if ( arr->nd != rank)
     throw copy_exception ()
      <<"   A deep copy of the object would be necessary while views are supposed to guarantee to present a *view* of the python data.\n"
      <<"   Rank mismatch . numpy array is of rank "<< arr->nd << "while you ask for rank "<< rank<<". \n";

    if ((Order == 'C') && (PyArray_ISFORTRAN(arr))) 
     throw copy_exception ()
      <<"   A deep copy of the object would be necessary while views are supposed to guarantee to present a *view* of the python data.\n"
      <<"     The numpy is in Fortran order while it is expected in C order. \n";

   if ((Order == 'F') && (!PyArray_ISFORTRAN(arr))) 
     throw copy_exception ()
      <<"   A deep copy of the object would be necessary while views are supposed to guarantee to present a *view* of the python data.\n"
      <<"     The numpy is not in Fortran order as it is expected. \n";

    numpy_obj = X; Py_INCREF(X); 
   }
   else { 
    bool ForceCast = false;// Unless FORCECAST is present in flags, this call will generate an error if the data type cannot be safely obtained from the object.
    int flags = (ForceCast ? NPY_FORCECAST : 0) ;// do NOT force a copy | (make_copy ?  NPY_ENSURECOPY : 0);
    if (!(PyArray_Check(X) && (Order=='D'))) flags |= (Order =='F' ? NPY_F_CONTIGUOUS : NPY_C_CONTIGUOUS); //impose mem order
    //if (!(PyArray_Check(X) && (Order=='D'))) flags |= (Order =='F' ? NPY_FARRAY : NPY_CARRAY); //impose mem order
    numpy_obj= PyArray_FromAny(X,PyArray_DescrFromType(elementsType), rank,rank, flags , NULL );

    // do several checks
    if (!numpy_obj) {// The convertion of X to a numpy has failed !
     if (PyErr_Occurred()) {PyErr_Print();PyErr_Clear();}
     TRIQS_RUNTIME_ERROR<<"numpy interface : the python object  is not convertible to a numpy. ";
    }
    assert (PyArray_Check(numpy_obj)); assert((numpy_obj->ob_refcnt==1) || ((numpy_obj ==X)));

    PyArrayObject *arr_obj = (PyArrayObject *)numpy_obj;
    try {
     if (arr_obj->nd!=rank)  TRIQS_RUNTIME_ERROR<<"numpy interface : internal error : dimensions do not match";
     if (arr_obj->descr->type_num != elementsType) 
      TRIQS_RUNTIME_ERROR<<"numpy interface : internal error : incorrect type of element :" <<arr_obj->descr->type_num <<" vs "<<elementsType;
     if (Order == 'F') { assert (PyArray_ISFORTRAN(numpy_obj));}
     else assert(PyArray_ISCONTIGUOUS(numpy_obj));
    }
    catch(...) { Py_DECREF(numpy_obj); throw;} // make sure that in case of problem, the reference counting of python is still ok...
   }

   PyArrayObject *arr_obj = (PyArrayObject *)numpy_obj;
   // extract strides and lengths
   const size_t dim =arr_obj->nd; 
   std::vector<size_t> lengths (dim);
   std::vector<std::ptrdiff_t> strides(dim); 
   for (size_t i=0; i< dim ; ++i) {
    lengths[i] = size_t(arr_obj-> dimensions[i]); 
    strides[i] = std::ptrdiff_t(arr_obj-> strides[i])/sizeof(ValueType);
   }
   // make the indexmap	
   IndexMapType ind (mini_vector<size_t,IndexMapType::rank>(lengths), mini_vector<std::ptrdiff_t,IndexMapType::rank>(strides),0);
   // make the storage (does not throw)
   storages::shared_block<ValueType> sto(numpy_obj,false);
   return std::make_pair(ind,sto);
  }

} // close namespace numpy_interface  

//-----------------------------------------------------

template <typename ValueType, int D, typename Opt> 
array<ValueType,D,Opt>::array (PyObject * X):impl_type(typename array<ValueType,D,Opt>::indexmap_type(),typename array<ValueType,D,Opt>::storage_type()) { 
 try { boost::tie(this->indexmap_,this->storage_) =  numpy_interface::extract_len_stri_block<indexmap_type,value_type>(X, true); }
 catch(numpy_interface::copy_exception ){ // intercept only this one...
  TRIQS_RUNTIME_ERROR<<"array construction from a numpy : the copy has failed for a mysterious reason ...";
 }
}

template <typename ValueType, typename Opt> 
matrix<ValueType,Opt>::matrix (PyObject * X):impl_type(typename matrix<ValueType,Opt>::indexmap_type(),typename matrix<ValueType,Opt>::storage_type()) { 
 try { boost::tie(this->indexmap_,this->storage_) =  numpy_interface::extract_len_stri_block<indexmap_type,value_type>(X, true); }
 catch(numpy_interface::copy_exception ){ // intercept only this one...
  TRIQS_RUNTIME_ERROR<<"matrix construction from a numpy : the copy has failed for a mysterious reason ...";
 }
}

template <typename ValueType, typename Opt> 
vector<ValueType,Opt>::vector (PyObject * X):impl_type(typename vector<ValueType,Opt>::indexmap_type(),typename vector<ValueType,Opt>::storage_type()) { 
 try { boost::tie(this->indexmap_,this->storage_) =  numpy_interface::extract_len_stri_block<indexmap_type,value_type>(X, true); }
 catch(numpy_interface::copy_exception ){ // intercept only this one...
  TRIQS_RUNTIME_ERROR<<"array construction from a numpy : the copy has failed for a mysterious reason ...";
 }
}


template <typename ValueType, int D, typename Opt> 
array_view<ValueType,D,Opt>::array_view (PyObject * X):impl_type(
  typename array_view<ValueType,D,Opt>::indexmap_type(),
  typename array_view<ValueType,D,Opt>::storage_type()) { 
 try { 
  boost::tie(this->indexmap_,this->storage_) = 
   numpy_interface::extract_len_stri_block<typename array_view<ValueType,D,Opt>::indexmap_type,typename array_view<ValueType,D,Opt>::value_type>(X, false);
 }
 catch(numpy_interface::copy_exception s){
  TRIQS_RUNTIME_ERROR<<"array_view construction from a numpy : could not take a view of type array_view <T,rank,Opt> with "
   <<"\n   T = "<< triqs::utility::typeid_name(ValueType())
   <<"\n   rank = "<< D
   <<"\n   Opt = "<< triqs::utility::typeid_name(Opt())
   <<"\nfrom the python object \n"<< numpy_interface::object_to_string(X) 
   <<"\nThe error was :\n "<<s.what();
 }
}

template <typename ValueType, typename Opt> 
matrix_view<ValueType,Opt>::matrix_view (PyObject * X):impl_type(
  typename matrix_view<ValueType,Opt>::indexmap_type(),
  typename matrix_view<ValueType,Opt>::storage_type()) { 
 try { 
  boost::tie(this->indexmap_,this->storage_) = 
   numpy_interface::extract_len_stri_block<typename matrix_view<ValueType,Opt>::indexmap_type,typename matrix_view<ValueType,Opt>::value_type>(X, false);
 }
 catch(numpy_interface::copy_exception s){
  TRIQS_RUNTIME_ERROR<<"matrix_view construction from a numpy : could not take a view of type matrix_view <T,Opt> with "
   <<"\n   T = "<< triqs::utility::typeid_name(ValueType())
   <<"\n   Opt = "<< triqs::utility::typeid_name(Opt())
   <<"\nfrom the python object \n"<< numpy_interface::object_to_string(X) 
   <<"\nThe error was :\n "<<s.what();
 }
}

template <typename ValueType, typename Opt> 
vector_view<ValueType,Opt>::vector_view (PyObject * X):impl_type(
  typename vector_view<ValueType,Opt>::indexmap_type(),
  typename vector_view<ValueType,Opt>::storage_type()) { 
 try { 
  boost::tie(this->indexmap_,this->storage_) = 
   numpy_interface::extract_len_stri_block<typename vector_view<ValueType,Opt>::indexmap_type,typename vector_view<ValueType,Opt>::value_type>(X, false);
 }
 catch(numpy_interface::copy_exception s){
  TRIQS_RUNTIME_ERROR<<"vector_view construction from a numpy : could not take a view of type vector_view <T,Opt> with "
   <<"\n   T = "<< triqs::utility::typeid_name(ValueType())
   <<"\n   Opt = "<< triqs::utility::typeid_name(Opt())
   <<"\nfrom the python object \n"<< numpy_interface::object_to_string(X) 
   <<"\nThe error was :\n "<<s.what();
 }
}

//-----------------------------------------------------

namespace numpy_interface  {

 template<typename T,int rank, typename Opt >
  PyObject * from_array_view(array_view<T,rank,Opt > const & A, bool copy=false) {
   _import_array();
   typedef typename array_view<T,rank,Opt >::value_type value_type;
   const int elementsType (numpy_to_C_type<typename boost::remove_const<value_type>::type>::arraytype);
   npy_intp dims[rank],  strides[rank];
   for(size_t i =0; i<rank; ++i) { dims[i] = A.indexmap().lengths()[i]; strides[i] = A.indexmap().strides()[i]*sizeof(value_type); }
   const value_type * data = A.data_start();
   //int flags = NPY_ARRAY_BEHAVED & ~NPY_ARRAY_OWNDATA;;// for numpy2
   int flags = NPY_BEHAVED & ~NPY_OWNDATA;
   PyObject* res  = PyArray_NewFromDescr(&PyArray_Type, PyArray_DescrFromType(elementsType), (int) rank, dims, strides, (void*) data,  flags, NULL);

   if (!res) { 
    if (PyErr_Occurred()) {PyErr_Print();PyErr_Clear();}
    TRIQS_RUNTIME_ERROR<<" array_view_from_numpy : the python numpy object could not be build";
   }
   if (!PyArray_Check(res)) TRIQS_RUNTIME_ERROR<<" array_view_from_numpy : internal error : the python object is not a numpy";
   PyArrayObject * arr = (PyArrayObject *)(res);
   //PyArray_SetBaseObject(arr,  A.storage().new_ref_to_guard());
   arr->base =  A.storage().new_ref_to_guard();

#ifdef TRIQS_ARRAYS_NUMPY_INTERFACE_DEBUG
   for(size_t i =0; i<rank; ++i) { cerr<<"i = " << i << PyArray_DIMS(arr)[i] << PyArray_STRIDES(arr)[i]<<endl;}
   cerr<<"nd = "<<arr->nd<<endl;
   cerr<<"flags = "<<arr->flags<<endl;
   cerr<< "data = "<< (void *)arr->data <<endl;
   cerr <<  ( (value_type*) arr->data)  [0] <<endl;
   cerr <<  ( (value_type*) arr->data)  [1] <<endl;
#endif

   assert( arr->flags == (arr->flags & ~NPY_OWNDATA));
   if (copy)  { 
    PyObject * na = PyObject_CallMethod(res,(char*)"copy",NULL);
    Py_DECREF(res);
    assert(((PyArrayObject *)na)->base ==NULL);
    res = na;
   }
   return res;
  }
 // array : make a deepcopy
 template<typename T,int rank, typename Opt >
  PyObject * from_array(array<T,rank,Opt > const & A) { return from_array_view ( array_view<T,rank,Opt > (A), true);}

}}}

#include <boost/python/handle.hpp>
#include <boost/python/slice.hpp>
#include "../../python_tools/converters.hpp"

// converter funtions
namespace triqs { namespace python_tools { 
 namespace Py_to_C { 

  template<>
   struct convert<triqs::arrays::range> { 
    static bool possible (boost::python::object obj) { return PySlice_Check(obj.ptr());} 
    static triqs::arrays::range invoke(boost::python::object obj) { 
     using boost::python::extract;
     PySliceObject* sl = (PySliceObject*)(obj.ptr());     
     return triqs::arrays::range(extract<long>(sl->start),extract<long>(sl->stop),extract<long>(sl->step));
    }
   };

  template<typename T,int rank, typename Opt >  
   struct convert<triqs::arrays::array<T,rank,Opt>  > { 
    static bool possible (boost::python::object obj) { 
     try  { triqs::arrays::array<T,rank,Opt>  A(obj.ptr()); }
     catch(triqs::runtime_error s) { std::cerr<<s.what()<<std::endl; return false;}
     return true;
    }
    static triqs::arrays::array<T,rank,Opt>  invoke(boost::python::object X) { return  triqs::arrays::array<T,rank,Opt> (X.ptr());}
   };

  template<typename T,int rank, typename Opt >  
   struct convert<triqs::arrays::array_view<T,rank,Opt>  > { 
    static bool possible (boost::python::object obj) { 
     try  { triqs::arrays::array_view<T,rank,Opt>  A(obj.ptr()); }
     catch(triqs::runtime_error s) { std::cerr<<s.what()<<std::endl; return false;}
     return true;
    }
    static triqs::arrays::array_view<T,rank,Opt>  invoke(boost::python::object X) { return  triqs::arrays::array_view<T,rank,Opt> (X.ptr());}
   };

#define CONVERT_MV(OO)\
  template<typename T, typename Opt >\
  struct convert<triqs::arrays::OO<T,Opt>  > {\
   static bool possible (boost::python::object obj) {\
    try  { triqs::arrays::OO<T,Opt>  A(obj.ptr()); }\
    catch(triqs::runtime_error s) { std::cerr<<s.what()<<std::endl; return false;}\
    return true;\
   }\
   static triqs::arrays::OO<T,Opt>  invoke(boost::python::object X) { return  triqs::arrays::OO<T,Opt> (X.ptr());}\
  };

  CONVERT_MV(matrix); 
  CONVERT_MV(matrix_view);
  CONVERT_MV(vector);
  CONVERT_MV(vector_view); 
#undef CONVERT_MV
 }

 namespace C_to_Py { 

  template<>
   struct convert<triqs::arrays::range> { 
    static boost::python::object invoke (triqs::arrays::range const & a) { return boost::python::slice(a.first(), a.last(), a.step());}
   };

  template<typename T,int rank, typename Opt >  
   struct convert<triqs::arrays::array<T,rank,Opt> > { 
    static boost::python::object invoke (triqs::arrays::array<T,rank,Opt>  const & a) { 
     PyObject * p = triqs::arrays::numpy_interface::from_array(a);
     return boost::python::object( boost::python::handle<>(p)); // p is a NEW reference ...
    }
   };

  template<typename T,int rank, typename Opt >  
   struct convert<triqs::arrays::array_view<T,rank,Opt> > { 
    static boost::python::object invoke (triqs::arrays::array_view<T,rank,Opt>  const & a) { 
     PyObject * p = triqs::arrays::numpy_interface::from_array_view(a);
     return boost::python::object( boost::python::handle<>(p)); // p is a NEW reference ...
    }
   };

  template<typename T, typename Opt >  
   struct convert<triqs::arrays::matrix<T,Opt> > { 
    static boost::python::object invoke (triqs::arrays::matrix<T,Opt> a) { //makes the copy here ...
     return convert<triqs::arrays::array_view<T,2,Opt> >::invoke(a);
    } 
   };

  template<typename T, typename Opt >  
   struct convert<triqs::arrays::matrix_view<T,Opt> > { 
    static boost::python::object invoke (triqs::arrays::matrix_view<T,Opt> a) { 
     return convert<triqs::arrays::array_view<T,2,Opt> >::invoke(a); 
    } 
   };

  template<typename T, typename Opt >  
   struct convert<triqs::arrays::vector<T,Opt> > { 
    static boost::python::object invoke (triqs::arrays::vector<T,Opt> a) { //makes the copy here ...
     return convert<triqs::arrays::array_view<T,1,Opt> >::invoke(a); 
    } 
   };

  template<typename T, typename Opt >  
   struct convert<triqs::arrays::vector_view<T,Opt> > { 
    static boost::python::object invoke (triqs::arrays::vector_view<T,Opt> a) { 
     return convert<triqs::arrays::array_view<T,1,Opt> >::invoke(a); 
    } 
   };
 }
}}
#endif
#endif
