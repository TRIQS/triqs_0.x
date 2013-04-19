
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

#ifndef PYTHON_ITERATORS_H
#define PYTHON_ITERATORS_H

#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <list>

namespace triqs { namespace python_tools { 

/** @defgroup Iteratorgroup python-C++ list/dict transcription tools.

    A few personnal iterators for python-C++ interface.

    boost python already contains necessary tools like stl_input_iterator.
    The following classes are just a bit more easy to use with a lighter
    syntax, but they are more specific.
*/
 /**
    @ingroup Iteratorgroup
    Iterates on a python list and transcribe the objects in C++
    
    Usage : 
     If L is a python::list of objects to be transcribed as T

     \code
     for (IteratorOnPythonList<T> p(L); !p.atEnd(); ++p)
       *p is a const & ref to the object T
     \endcode
   */
  template<typename T>
  class IteratorOnPythonList { 
    boost::python::stl_input_iterator<boost::python::object> beg,end;
    std::list<boost::python::object> l;
    std::list<boost::python::object>::const_iterator p;
    T result;
    inline void convert() {result = boost::python::extract<T>((*p));}    
  public:
    ///
    IteratorOnPythonList( const boost::python::list &L): beg(L),end(),l(beg,end) { p=l.begin(); }
    ///
    IteratorOnPythonList& operator++() {++p;  return *this;}
    ///
    const T & operator*(){ convert(); return (result);}
    ///
    bool atEnd () const { return (p==l.end());} 
  };
 
  /**
    @ingroup Iteratorgroup
    Iterates on a python list of 2-tuples and transcribe the objects of the tuples into T1 and T2 C++ objects 
    respectively.
    
    Usage : 
     If A is a python list of (x1,x2)
     which are to be transcribed into tuple (T1,T2)

     \code
     for (IteratorOnPythonListOf2Tuples<T1,T2> p(L); !p.atEnd(); ++p) {
       p->x1 is a const & on x1 
       p->x2 is a const & on x2
     \endcode
    
   */
  template<typename T1, typename T2>
  class IteratorOnPythonListOf2Tuples { 
  public : 
    struct pair {
      /// first element
      T1 x1; 
      /// second element
      T2 x2;};
  private:
    boost::python::stl_input_iterator<boost::python::object> beg,end;
    std::list<boost::python::object> l;
    std::list<boost::python::object>::const_iterator p;
    pair result;
    inline void convert() {result.x1 = boost::python::extract<T1>((*p)[0]);result.x2 = boost::python::extract<T2>((*p)[1]);}
  public:
    ///
    IteratorOnPythonListOf2Tuples( const boost::python::list &L): beg(L),end(),l(beg,end){ p=l.begin(); }
    ///
    IteratorOnPythonListOf2Tuples& operator++() {++p;  return *this;}
    ///
    const pair * operator->(){ convert(); return &(result);}
    ///
    const pair & operator*() { convert(); return (result);}
    /// Is the iterator at the end ?
    bool atEnd () const { return (p==l.end());} 
  };
  
   /**
     @ingroup Iteratorgroup
    Iterates on a python list of 3-tuples and transcribe the objects of the tuples into T1, T2 and T3 C++ objects 
    respectively.
    
    Usage : 
     If A is a python list of (x1,x2,x3)
     which are to be transcribed into tuple (T1,T2,T3)

     \code
     for (IteratorOnPythonListOf3Tuples<T1,T2,T3> p(L); !p.atEnd(); ++p) {
       p->x1 is a const & on x1 
       p->x2 is a const & on x2
       p->x3 is a const & on x3
     \endcode
    
   */
  template<typename T1, typename T2, typename T3>
  class IteratorOnPythonListOf3Tuples { 
  public : 
    struct pair {
      /// first element
      T1 x1; 
      /// second element
      T2 x2; 
      /// third element
      T3 x3;};
  private:
    boost::python::stl_input_iterator<boost::python::object> beg,end;
    std::list<boost::python::object> l;
    std::list<boost::python::object>::const_iterator p;
    pair result;
    inline void convert() {
      result.x1 = boost::python::extract<T1>((*p)[0]);
      result.x2 = boost::python::extract<T2>((*p)[1]);
      result.x3 = boost::python::extract<T3>((*p)[2]);}
  public:
    ///
    IteratorOnPythonListOf3Tuples( const boost::python::list &L): beg(L),end(),l(beg,end){ p=l.begin(); }
    ///
    IteratorOnPythonListOf3Tuples& operator++() {++p;  return *this;}
    ///
    const pair * operator->(){ convert(); return &(result);}
    ///
    const pair & operator*() { convert(); return (result);}
    /// Is the iterator at the end ?
    bool atEnd () const { return (p==l.end());} 
  };
  

  /**
    @ingroup Iteratorgroup
    Iterates on a python dict and transcribe the key and values into T1 and T2 C++ objects 
    respectively.
    
     Usage : 
     If A is a python::dict of (key,values)
     which are to be transcribed into tuple (T1,T2)

     \code
     for (IteratorOnPythonDict<T1,T2> p(L); !p.atEnd(); ++p)
       p->key is a const & on the key 
       p->val is a const & on the value 
     \endcode
    
   */
  template<typename T1, typename T2>
  class IteratorOnPythonDict { 
  public : 
    /// To store the result
    struct pair { 
      /// Key 
      T1 key; 
      /// Value
      T2 val;};
  private:
    boost::python::stl_input_iterator<boost::python::object> beg,end;
    std::list<boost::python::object> l;
    std::list<boost::python::object>::const_iterator p;
    pair result;
    inline void convert() {result.key = boost::python::extract<T1>((*p)[0]);result.val = boost::python::extract<T2>((*p)[1]);}
  public:
    ///
    IteratorOnPythonDict( const boost::python::dict &d):beg(d.items()),end(),l(beg,end) { p=l.begin(); }
    ///
    IteratorOnPythonDict& operator++() {++p;  return *this;}
    ///
    const pair * operator->(){ convert(); return &(result);}
    ///
    const pair & operator*(){ convert(); return (result);}
    /// Is the iterator at the end ?
    bool atEnd () const { return (p==l.end());} 
  };
  
} }

#endif
