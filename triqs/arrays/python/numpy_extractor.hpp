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
#ifndef TRIQS_ARRAYS_NUMPY_EXTRACTOR_H
#define TRIQS_ARRAYS_NUMPY_EXTRACTOR_H 

#ifdef TRIQS_WITH_PYTHON_SUPPORT
#include "../storages/shared_block.hpp"
#include "triqs/utility/exceptions.hpp"
#include "numpy/arrayobject.h"

namespace triqs { namespace arrays { namespace numpy_interface  {

 //inline std::string object_to_string (const bpy::object & O1) { return bpy::extract<std::string>(bpy::str(O1)); }
 //inline std::string object_to_string (PyObject * p) { bpy::object obj ( bpy::borrowed (p)); return object_to_string(obj); }
 
 inline std::string object_to_string (PyObject * p) { 
  if (!PyString_Check(p)) TRIQS_RUNTIME_ERROR<<" Internal error, expected a python string .....";
  return PyString_AsString(p); 
 } 

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

 template<typename IndexMapType, typename ValueType > struct numpy_extractor { 

  numpy_extractor ( PyObject * X, bool allow_copy) {
   if (X==NULL) TRIQS_RUNTIME_ERROR<<"numpy interface : the python object is NULL !";
   if (_import_array()!=0) TRIQS_RUNTIME_ERROR <<"Internal Error in importing numpy";

   // make sure IndexMap is cuboid ?s
   static const char Order = IndexMapType::index_order_type::C_or_F;
   static_assert( ((Order=='F')||(Order=='C')), "Ordering must be C or Fortran");

   if ((!PyArray_Check(X)) && (Order=='D'))
    TRIQS_RUNTIME_ERROR<<"numpy interface : the python object is not a numpy and you ask me to deduce the ordering in memory !";

   const int elementsType (numpy_to_C_type<typename boost::remove_const<ValueType>::type>::arraytype);
   int rank = IndexMapType::rank;
   static const char * error_msg = "   A deep copy of the object would be necessary while views are supposed to guarantee to present a *view* of the python data.\n";

   if (!allow_copy) {
    // in case of a view, we decide ourselves if we can do it...
    // a previous uses PyArray_FromAny, but behaviour changes between different version of numpy, so it is not portable...
    if (!PyArray_Check(X)) 
     throw copy_exception () << error_msg<<"   Indeed the object was not even an array !\n";  

    if ( elementsType != PyArray_TYPE((PyArrayObject*)X)) 
     throw copy_exception () << error_msg<<"   The deep copy is caused by a type mismatch of the elements. \n";

    PyArrayObject *arr = (PyArrayObject *)X;    
    if ( arr->nd != rank)
     throw copy_exception () << error_msg<<"   Rank mismatch . numpy array is of rank "<< arr->nd << "while you ask for rank "<< rank<<". \n";

    if ((Order == 'C') && (PyArray_ISFORTRAN(arr))) 
     throw copy_exception () << error_msg<<"     The numpy is in Fortran order while it is expected in C order. \n";

    if ((Order == 'F') && (!PyArray_ISFORTRAN(arr))) 
     throw copy_exception () << error_msg<<"     The numpy is not in Fortran order as it is expected. \n";

    numpy_obj = X; Py_INCREF(X); 
   }
   else { 
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
     if (Order == 'F') { if (!PyArray_ISFORTRAN(numpy_obj)) TRIQS_RUNTIME_ERROR<<"numpy interface : internal error : should be Fortran array";}
     else {if (!PyArray_ISCONTIGUOUS(numpy_obj)) TRIQS_RUNTIME_ERROR<<"numpy interface : internal error : should be contiguous";}
    }
    catch(...) { Py_DECREF(numpy_obj); throw;} // make sure that in case of problem, the reference counting of python is still ok...
   }

   arr_obj = (PyArrayObject *)numpy_obj;
  } 

  ///
  IndexMapType indexmap() const {
   // extract strides and lengths
   const size_t dim =arr_obj->nd; 
   std::vector<size_t> lengths (dim);
   std::vector<std::ptrdiff_t> strides(dim); 
   for (size_t i=0; i< dim ; ++i) {
    lengths[i] = size_t(arr_obj-> dimensions[i]); 
    strides[i] = std::ptrdiff_t(arr_obj-> strides[i])/sizeof(ValueType);
   }
   return IndexMapType (mini_vector<size_t,IndexMapType::rank>(lengths), mini_vector<std::ptrdiff_t,IndexMapType::rank>(strides),0);
  }

  ///
  storages::shared_block<ValueType> storage() const { return storages::shared_block<ValueType> (numpy_obj,false); }

  private:
  PyObject * numpy_obj;
  PyArrayObject *arr_obj;
 };
}}} 
#endif
#endif
