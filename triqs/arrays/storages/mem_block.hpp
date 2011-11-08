
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

#ifndef TRIQS_MEM_BLOCK_H 
#define TRIQS_MEM_BLOCK_H
#include "../../utility/exceptions.hpp"
#include <boost/noncopyable.hpp>
#include <Python.h>
#include <numpy/arrayobject.h>
#include "./common.hpp"
#ifdef TRIQS_ARRAYS_DEBUG_TRACE_MEM
#include <iostream>
#endif

namespace triqs { namespace arrays { namespace storages { namespace details { 

 template<typename ValueType >
  class mem_block : boost::noncopyable {

   typedef typename boost::remove_const<ValueType>::type non_const_value_type;
   size_t size_;
   non_const_value_type * restrict p;
   PyObject * py_obj;
   static void delete_pointeur( void *ptr ) { 
#ifdef TRIQS_ARRAYS_DEBUG_TRACE_MEM
    std::cerr<< "deleting data block"<<(non_const_value_type*) ptr<< std::endl; 
#endif
    delete [] ( (non_const_value_type*) ptr) ;
   }   

   static void import_numpy_array() { 
    static bool init = false;
    //if (init) { std::cerr<<" memblock import arry : already done"<<std::endl; return;}
    int _r = _import_array();assert(_r==0);
    //std::cerr<<" memblock import array done"<<std::endl;
    init = true;
   }

   public : 

   mem_block():size_(0),p(NULL),py_obj(NULL){}

   mem_block (size_t s):size_(s),p(new non_const_value_type[s]),py_obj(NULL){
#ifdef TRIQS_ARRAYS_DEBUG_TRACE_MEM
    std::cerr<< " construct meme block from C :p ="<<p<< "  py_obj = "<< py_obj<<std::endl;
#endif
   }

   mem_block (PyObject * obj, bool borrowed ) { 
#ifdef TRIQS_ARRAYS_DEBUG_TRACE_MEM
    std::cerr<< " construct meme block from pyobject"<<obj<< " # ref ="<<obj->ob_refcnt<<" borrowed = "<< borrowed<<std::endl; 
#endif 
    mem_block::import_numpy_array(); 
    //std::cerr<<"  ref o obj = "<< obj->ob_refcnt<<std::endl;
    if (borrowed) Py_INCREF(obj);
    py_obj = obj;
    //assert(obj);
    if ( (!obj) || (!PyArray_Check(obj))) TRIQS_RUNTIME_ERROR<<"Internal error : mem_block construct from pyo : obj is not an array";
    PyArrayObject * arr = (PyArrayObject *)(obj);
    size_ = PyArray_SIZE(arr);
    this->p = (non_const_value_type*)PyArray_DATA(arr); 
   }

   void activate_python() { 
#ifdef TRIQS_ARRAYS_DEBUG_TRACE_MEM
    std::cerr<< " activating python "<<p<< "  py_obj = "<< py_obj<< std::endl;
#endif
    if (py_obj==NULL) { py_obj = PyCObject_FromVoidPtr( (void*) p, &mem_block<ValueType>::delete_pointeur);} 
   }

   // delete memory manually iif py_obj is not set. Otherwise the python interpreter will do that for us.
   ~mem_block(){ 
#ifdef TRIQS_ARRAYS_DEBUG_TRACE_MEM
    std::cerr<< "deleting mem block p ="<<p<< "  py_obj = "<< py_obj<< std::endl;
    if (py_obj) std::cerr<<"    ref of py obj"<<py_obj->ob_refcnt<<std::endl;
#endif
    if (py_obj==NULL) {
     if (p) {
#ifdef TRIQS_ARRAYS_DEBUG_TRACE_MEM
      std::cerr<< "deleting data block 2"<<std::endl; 
#endif
      delete[] p;
     }
    } 
    else Py_DECREF(py_obj); 
   } 

   void operator=(const mem_block & X) {
    assert( py_obj==NULL); assert(size_==X.size_);assert(p); assert(X.p);
    memcpy (p,X.p,size_ * sizeof(ValueType));
   }

   mem_block * clone () const  { 
#ifdef TRIQS_ARRAYS_DEBUG_TRACE_MEM
    std::cerr<< "cloning "<<std::endl; 
#endif
    mem_block * r= new mem_block(size_); 
    if ((py_obj==NULL) || (PyCObject_Check(py_obj))) { (*r) = (*this); return r;}
    // else make a new copy of the numpy ...
    mem_block::import_numpy_array(); 
    assert(PyArray_Check(py_obj));
    //std::cerr<< " python cloning "<<std::endl; 
    if ( ( PyArray_ISFORTRAN(py_obj)) || (PyArray_ISCONTIGUOUS(py_obj)))  { 
     memcpy (r->p,PyArray_DATA(py_obj),size_ * sizeof(ValueType));
    }
    else { // if the py_obj is not contiguous, first let numpy copy it properly
     PyObject * na = PyObject_CallMethod(py_obj,(char *)"copy",NULL);
     assert(na); assert( ( PyArray_ISFORTRAN(na)) || (PyArray_ISCONTIGUOUS(na)));
     memcpy (r->p,PyArray_DATA(na),size_ * sizeof(ValueType));
     Py_DECREF(na);
    }
    return r;
   }

   PyObject * new_ref_to_guard() { activate_python(); Py_INCREF(py_obj); return py_obj;} 

   non_const_value_type & operator[](size_t i) {return p[i];}

   size_t size() const {return size_;}

   template<class Archive>
    void save(Archive & ar, const unsigned int version) const { 
     ar << boost::serialization::make_nvp("size",size_);
     for (size_t i=0; i<size_; ++i) ar << boost::serialization::make_nvp("data",p[i]); 
    }

   template<class Archive>
    void load(Archive & ar, const unsigned int version) { 
     ar >> size_;
     assert (p==NULL); 
     p = new non_const_value_type[size_]; 
     for (size_t i=0; i<size_; ++i) ar >> p[i]; 
    }
   BOOST_SERIALIZATION_SPLIT_MEMBER();
  };
}}}}//namespace triqs::arrays 
#endif

