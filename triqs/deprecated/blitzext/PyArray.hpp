
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

#ifndef PYARRAY_237849273984_H
#define PYARRAY_237849273984_H

#include <Python.h>
#include <triqs/utility/report_stream.hpp>
#include <numpy/arrayobject.h>
#include <blitz/array.h>
#include <boost/python.hpp>
#include <string>
#include <complex>

 #define FATAL(s) { std::stringstream fs_for_fatal; fs_for_fatal<<s; throw fs_for_fatal.str();}

/** @cond INTERNAL */
namespace PyArray_internal {
  
  // Transcription of numpy types to C++ types
  template <class T> struct PyArrayType {};
  template <> struct PyArrayType<int> {static const int  arraytype = NPY_INT;};
  template <> struct PyArrayType<long> {static const int arraytype = NPY_LONG;};
  template <> struct PyArrayType<double> {static const int arraytype = NPY_DOUBLE;};
  template <> struct PyArrayType<std::complex<double> > {static const int arraytype = NPY_CDOUBLE;};
  template <> struct PyArrayType<std::string > {static const int arraytype = NPY_STRING;};
    
  // default value for various types
  template <class T> inline  T TheDefaults() {return T();}
  template <> inline bool TheDefaults<bool>() {return false;}
  template <> inline int TheDefaults<int>() {return 0;}
  template <> inline long TheDefaults<long>() {return 0;}
  template <> inline double TheDefaults<double>() {return 0;}
  template <> inline std::complex<double> TheDefaults<std::complex<double> >() {return 0;}
  template <> inline std::string TheDefaults<std::string>() {return "";}
}
/** @endcond  */

/// Specify the ordering possibilities of an array
enum PyArrayOrderingInMemory {FortranOrder, COrder,SameOrder};
//enum PyArrayOrderingInMemoryFC {Fortran, C};

/// CopyPolicy
enum PyArrayCopyPolicy {ViewIfPossible, CopyMandatory, AssertView}; 

/** @cond INTERNAL */
namespace PyArray_internal {

  // This is a technical class that construct the numpy part of the array. Internal Use only.
  ///
  template<class T,int N> class PyArray_ancestor {
  public:
    PyArray_ancestor (PyObject *x,std::string mess, 
		      PyArrayCopyPolicy copyPolicy, PyArrayOrderingInMemory Order, 
		      bool ForceCast , bool AdjustDimensionNumber ):
      myshape(0),strides(0)
    { 
      int _r = _import_array();assert(_r==0);
      // If x is an array : 
      //   - if Order is same, don't change it
      //   - else impose it (may provoque a copy).
      // if x is not array : Order = FortranOrder or SameOrder - > Fortran order otherwise C
      if (PyArray_Check(x) && (Order==SameOrder))
	flags = (ForceCast ? NPY_FORCECAST : 0) | (copyPolicy==CopyMandatory ?  NPY_ENSURECOPY : 0);
      else
	flags = (Order == FortranOrder ? NPY_FARRAY : NPY_CARRAY) | (ForceCast ? NPY_FORCECAST : 0) | (copyPolicy==CopyMandatory ?  NPY_ENSURECOPY : 0);
      // cout<<" REFS CO|NT AVANT "<<x->ob_refcnt<<"  P="<<x<<endl;
       arr_py= (AdjustDimensionNumber ? 
	       PyArray_FromAny(x,PyArray_DescrFromType(PyArray_internal::PyArrayType<T>::arraytype),N,N, flags , NULL ) : 
	       PyArray_FROM_OTF(x,PyArray_internal::PyArrayType<T>::arraytype, flags ) );
      //cout<<" REFS CO|NT APRES "<<arr_py->ob_refcnt<<"  P="<<arr_py<<endl;
      if (!arr_py) {
	if (PyErr_Occurred()) {PyErr_Print();PyErr_Clear();}
	//FATAL("PyArray <Array<T,n>  : the python object *"<<string(extract<string>(str(object(handle<>(arr_py)))))<<" is not convertible to a numpy ... "<<mess);
	FATAL("PyArray <Array<T,n>  : the python object is not convertible to a numpy ... "<<mess);
      }

      //cout<<"ORDER"<<Order<<endl;
      if (Order == FortranOrder) { //cout<<"MADE FORTRAN"<<endl; 
	assert (PyArray_ISFORTRAN(arr_py));}

      //      cout<<"COPY policy"<< copyPolicy<<  "  "<<(copyPolicy==AssertView)<<endl;

      if ( (copyPolicy==AssertView) && (arr_py !=x)) 
	FATAL("Assert View not respected in PyArray !");
      setup(mess);


#ifdef PYARRAY_DEBUG
      if (arr_py !=x) std::cout<<"COPY MADE"<<std::endl;
#endif
    }
    
    PyArray_ancestor (const blitz::TinyVector<npy_intp,N> & dim, bool fortran ): myshape(0),strides(0)
    { 
      int _r = _import_array();assert(_r==0);
      npy_intp dims[N]; for (int i=0;i<N;i++) dims[i] = dim[i];

      // PyObject *res = PyArray SimpleNew (N, dims, PyArrayType<T>::arraytype);
      arr_py = (fortran ? 
		PyArray_New(&PyArray_Type, N, dims, PyArray_internal::PyArrayType<T>::arraytype, NULL, NULL, 0, NPY_FARRAY ,NULL) : 
		PyArray_New(&PyArray_Type, N, dims, PyArray_internal::PyArrayType<T>::arraytype, NULL, NULL, 0, 0 ,NULL));

      if (!arr_py) {
	if (PyErr_Occurred()) {PyErr_Print();PyErr_Clear();}
	FATAL("PyArray <Array<T,n>  : can not construct ");
      }

      PyArrayObject *arr_obj = (PyArrayObject *)arr_py;
      if (arr_obj->descr->type_num != PyArrayType<T>::arraytype) 
	std::cout<<" CE N'EST PAS LE TYPE DEMANDE "<<  arr_obj->descr->type_num << "   "<< PyArrayType<T>::arraytype<<std::endl;
      // Bug Bizarre : numpy 1.0.3 NPY_CARRAY -> give me a F_Array !
      setup(" INTERNAL PB setup");
    };

 
    PyArray_ancestor( const PyArray_ancestor<T,N> & x,bool copy) : myshape(x.myshape),strides(x.strides),
								   flags(x.flags),arr_py(x.arr_py)
    {
      int _r = _import_array();assert(_r==0);
      flags = (copy ?  NPY_FORCECAST | NPY_ENSURECOPY : 0);
      if (!PyArray_Check(arr_py)) FATAL("PyArray <Array<T,n> : copy constructor : First arg is not a Numeric array ");
      arr_py= PyArray_FROM_OTF(arr_py,PyArray_internal::PyArrayType<T>::arraytype, flags);
	   
      if (!arr_py) {
	std::cout<<" AN ERORR HAS OCCURED"<<std::endl;
	if (PyErr_Occurred()) {PyErr_Print();PyErr_Clear();}
	//FATAL("PyArray <Array<T,n>  : the python object **"<<string(extract<string>(str(object(handle<>(arr_py)))))<<" is not convertible to a numpy ... ");
	FATAL("PyArray <Array<T,n> : copy constructor  : the python object provided is not convertible to a numpy ... ");
      }
      if ((copy==false) && (arr_py !=x.arr_py)) FATAL("   PyArray_ancestor : copy NOT requested but made by numpy !!");
      setup("");
    }

    ~PyArray_ancestor(){
      if (arr_py)  {
#ifdef PYARRAY_DEBUG
	std::cout<<" REFS AVANT DESTRCUTION "<<arr_py->ob_refcnt<<"  P="<<arr_py<<std::endl;
#endif
	Py_DECREF(arr_py);} 
    }
    
  protected : 
    blitz::TinyVector<int,N>  myshape,strides;
    int flags;
    PyObject* arr_py;
    
  private : 
    void setup(std::string mess) {	
      int _r = _import_array();assert(_r==0);
      PyObject * obj = arr_py;
      if (!PyArray_Check(obj)) FATAL("Blitz_View_of_Numpy : First arg is not a Numeric array "<<mess);
      PyArrayObject *arr_obj = (PyArrayObject *)obj;

#ifdef PYARRAY_DEBUG
      if (arr_obj->descr->type_num != PyArrayType<T>::arraytype) {
	std::cout<<"PBPB "<< arr_obj->descr->type_num << "  "<< PyArrayType<T>::arraytype<<std::endl;
	std::cout<< PyArrayType<int>::arraytype<< "   "<< PyArrayType<long>::arraytype<<std::endl;
      }
#endif

      if (arr_obj->descr->type_num != PyArrayType<T>::arraytype) 
	FATAL("PyArray Construction : The Numeric array does not contain the correct type of number"<<mess);
      if (arr_obj->nd!=N) FATAL("PyArray Construction :  dimensions do not match"<<mess);
      int T_size = sizeof(T);
      npy_intp *arr_dimensions = arr_obj-> dimensions, *arr_strides = arr_obj-> strides;
      for (int i=0;i<N;++i) {
	myshape[i]   = int(arr_dimensions[i]);
	strides[i] = int(arr_strides[i])/T_size;
      }
    }
  };
} //internal namespace
/** @endcond  */


  /**
     Class PyArray : view of numpy array as C++ blitz++ array.
     Given a PyObject to a numpy, provides a VIEW in this numpy array as a blitz array if possible, 
     a view on the copy when it is not possible. In any case, a view on a python numpy.
     The promotion of the type is automatic when possible.
  */
template<class T,int N> 
class PyArray : protected PyArray_internal::PyArray_ancestor<T,N>, public  blitz::Array<T,N> { 
public:
  /**
     Takes the PyObject, and makes a numpy.array out of it of the type requested.
     Equivalent to doing numpy.array(x, copy=True, order = ...) in python.
     
     If possible, it will return a view on the array.
     If not possible, it will create a new numpy and return a view on it.
     Indeed, it is not always possible to make a simple, e.g. you ask PyArray<double,1> (X)
     where X is a numpy of int.
     
     Parameters are : 
       - CopyPolicy : 
           - ViewIfPossible -> return a view of the Pyobject if possible (see above).
           - CopyMandatory  -> force the copy of the python object before making the view
	   - AssertView : if a view is not possible, raise a exception.
       - Order : FortranOrder, COrder,SameOrder
           - SameOrder -> Keep the same order as in numpy 
	   - FortranOrder -> Force Fortran order (may require copy)
	   - COrder -> Force C order (may require copy)
       - ForceCast [expert only] : Force the cast of the elements of the array to T. 
          If false, cast only occur when safe (cf numpy book) otherwise an error is raised.
       - AdjustDimensionNumber : if true, allows to build the n dim object from n-x dims ones.
  */
  PyArray(PyObject *x,std::string mess="", PyArrayOrderingInMemory Order= SameOrder, 
	  PyArrayCopyPolicy copyPolicy = ViewIfPossible, bool ForceCast = false, bool AdjustDimensionNumber = false ) :
    PyArray_internal::PyArray_ancestor<T,N>(x,mess,copyPolicy,Order,ForceCast,AdjustDimensionNumber),
    blitz::Array<T,N>( (T*) (PyArray_DATA(this->arr_py)),this->myshape,this->strides,blitz::neverDeleteData,
		       ( PyArray_ISFORTRAN(this->arr_py) ?    blitz::FortranArray<N>() : blitz::GeneralArrayStorage<N>()))  
  {};

  /**
     Equivalent to PyArray(x.ptr(), .....)
  */
  explicit PyArray(boost::python::object x,std::string mess="", PyArrayOrderingInMemory Order= SameOrder, 
	  PyArrayCopyPolicy copyPolicy = ViewIfPossible, bool ForceCast = false, bool AdjustDimensionNumber = false ) :
    PyArray_internal::PyArray_ancestor<T,N>(x.ptr(),mess,copyPolicy,Order,ForceCast,AdjustDimensionNumber),
    blitz::Array<T,N>( (T*) (PyArray_DATA(this->arr_py)),this->myshape,this->strides,blitz::neverDeleteData,
		       ( PyArray_ISFORTRAN(this->arr_py) ?    blitz::FortranArray<N>() : blitz::GeneralArrayStorage<N>()))  
  {};
    
/** @cond INTERNAL */
  /**
     Construct a new numpy array of type T, and dimension N.
     The array is initialized to 0.
   */
  #define C(NN)\
   explicit PyArray(BOOST_PP_ENUM_PARAMS(NN, int n), PyArrayOrderingInMemory order ):\
    PyArray_internal::PyArray_ancestor<T,N>( blitz::TinyVector<npy_intp,N>( BOOST_PP_ENUM_PARAMS(NN, n)) ,(order==FortranOrder)),\
    blitz::Array<T,N>( (T*) (PyArray_DATA(this->arr_py)),this->myshape,this->strides,blitz::neverDeleteData,\
		       ( PyArray_ISFORTRAN(this->arr_py) ?    blitz::FortranArray<N>() : blitz::GeneralArrayStorage<N>()))\
  { assert (N==NN);*this =0; }

  C(1);
  C(2);
  C(3);
  C(4);
  C(5);
#undef C
/** @endcond  */

  // This will provide all constructor, but compilation will only be ok iif NN=N

  /**
     Returns a VIEW of the current array unless copy=true
   */
  PyArray (const PyArray<T,N> & X, bool copy=false): 
    PyArray_internal::PyArray_ancestor<T,N>(X,copy),
    blitz::Array<T,N>( (T*) (PyArray_DATA(this->arr_py)),this->myshape,this->strides,blitz::neverDeleteData,
		       ( PyArray_ISFORTRAN(this->arr_py) ?    blitz::FortranArray<N>() : blitz::GeneralArrayStorage<N>()))  
  {};

  /**
     Returns a COPY of the current array
  */
  PyArray<T,N> copy() const { return PyArray(*this,true);}

  /**
     Same as for Array<T,N> : copy the data of X into *this
   */
  inline PyArray & operator= (const PyArray<T,N> & X) { 
    (*this) = (blitz::Array<T,N>) X; 
    return *this;
  }

  /** Same as for usual Array*/
  template <class R>
  inline  PyArray & operator = (const R & r) { blitz::Array<T,N>::operator=(r); return *this;}
	 
  /**
     Returns *this as a python::object (no copy made, ever).
  */
  boost::python::object as_BoostObject() const
  { return boost::python::object ( boost::python::borrowed( this->arr_py));}

  /**
     Returns a NEW reference to the array numpy
   */
  PyObject * as_PyObjectPtrNewRef() const { Py_INCREF(this->arr_py); return this->arr_py;}

  /**
     Returns a BORROWED reference to the array numpy
   */
  PyObject * as_PyObjectPtrBorrowedRef() const { return this->arr_py;}

  /**
     Cast itself as a  blitz::Array<T,N> explicitely
  */
  const blitz::Array<T,N> & as_blitzArray() const {return blitz::Array<T,N>(*this);}

  /**
     Cast itself as a  blitz::Array<T,N> explicitely
  */
  blitz::Array<T,N> & as_blitzArray() {return blitz::Array<T,N>(*this);}

};


//  Specific case of array of string. Not implemented.
template<int n> 
class PyArray<std::string,n> :  public  blitz::Array<std::string,n> { 
private:
  PyArray(PyObject *x,std::string mess="") { assert (0);}
};


/*----------------------------------------
  BOOST CONVERTER for numpy <-> Blitz++ arrays
  ------------------------------------------*/
 
/** @cond INTERNAL */
namespace PyArray_internal{
  
  template<class T, int n>
  struct PyArray_to_python { 
    static PyObject* convert (const PyArray<T,n> & a) { return a.as_PyObjectPtrNewRef(); }  
  };

  template<class T, int n>
  struct PyArray_from_python_str
  {
    PyArray_from_python_str(){
      boost::python::converter::registry::push_back(&convertible, &construct, boost::python::type_id< PyArray<T,n> >());
    }

    static void* convertible(PyObject* obj_ptr) {
      //int _r = _import_array();assert(_r==0);
      // if (!PyArray_Check(obj_ptr)) return 0;
       try  {
	PyArray<T,n> tmp (obj_ptr);
      }
      catch(std::string s) { return 0;} 
      return obj_ptr;
    }

    static void construct( PyObject* obj_ptr,
			   boost::python::converter::rvalue_from_python_stage1_data* data) {
      typedef boost::python::converter::rvalue_from_python_storage<PyArray<T,n> > storage_t;
      storage_t* the_storage = reinterpret_cast<storage_t*>( data );
      void* memory_chunk = the_storage->storage.bytes;
      //PyArray<T,n> * A =;
      new (memory_chunk) PyArray<T,n>(obj_ptr,"from converter: internal pb",SameOrder,ViewIfPossible,false,false);
      data->convertible = memory_chunk;
    }
  };
}
/** @endcond  */

template<class T, int n>
void register_converter_python_array() {
  boost::python::to_python_converter< PyArray<T,n>, PyArray_internal::PyArray_to_python<T,n> >();
  PyArray_internal::PyArray_from_python_str<T,n>();
}



template <class T>
class myPythonIteratorOnPyArray { 
public : 
  myPythonIteratorOnPyArray (const PyArray<T,1> & a) :A(a.as_PyObjectPtrBorrowedRef()),i(A.lbound(0)), stop(a.ubound(0)+1){}
  myPythonIteratorOnPyArray __iter__(){ return *this;}
  T next(){
    if (i==stop) {
      PyErr_SetString(PyExc_StopIteration, "");
      boost::python::throw_error_already_set();
    }
    T r(A(i)); i++; return r;}
private :
  const PyArray<T,1>  A;
  int i,stop;

};

#undef FATAL
#endif
