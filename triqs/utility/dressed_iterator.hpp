/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013 by O. Parcollet
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
#ifndef TRIQS_UTILITY_ITERATOR_DRESSING_H
#define TRIQS_UTILITY_ITERATOR_DRESSING_H
#include "./first_include.hpp"
#include <boost/iterator/iterator_facade.hpp>

namespace triqs { namespace utility {

 /**
  * Usage :
  *
  * Suppose we have an iterator it of type it_type and a little struct 
  * that views the data pointed to by it using a set of references, 
  * like : 
  *
  * struct _dress  {double const & tau; int & a;  _dress ( it_type const & it); };
  *
  * then dressed_iterator<it_type, _dress> is :
  *
  *  - a bidirectional iterator, STL compliant (?).
  *  - that is dereferenced to a _dress.
  *
  *  Example : See doc or test/utility/iterator_dressing.cpp ...
  *
  *
  * Variant : one can keep in the iterator a pointer to anything (typically the parent container)
  * which is useful in some cases
  * See iterator_dressing3 for an example
  */
 template< typename BidirectionalIteratorType, typename Dressing, typename DressingAuxiliaryArgumentPtrType = void> struct dressed_iterator;

 // specialization when a aux data is present
 template< typename BidirectionalIteratorType, typename Dressing, typename DressingAuxiliaryArgumentPtrType> 
  struct dressed_iterator : public boost::iterator_facade<dressed_iterator<BidirectionalIteratorType,Dressing,DressingAuxiliaryArgumentPtrType>, 
  Dressing, boost::bidirectional_traversal_tag,Dressing > {
   public : 
    dressed_iterator () {}
    dressed_iterator (BidirectionalIteratorType const & it, DressingAuxiliaryArgumentPtrType * aux): _it(it), _aux(aux) {}
    BidirectionalIteratorType const & get() const { return _it;}
    BidirectionalIteratorType       & get()       { return _it;}
    operator BidirectionalIteratorType() const { return _it;}
    DressingAuxiliaryArgumentPtrType * get_aux() { return _aux;}
    const DressingAuxiliaryArgumentPtrType * get_aux() const { return _aux;}
   private:
    friend class boost::iterator_core_access;
    void increment(){ ++_it;}
    void decrement(){ --_it;}
    bool equal(dressed_iterator const & other) const { return(other._it==_it);}
    Dressing dereference() const { return Dressing(_it,_aux); } 
    BidirectionalIteratorType _it;
    DressingAuxiliaryArgumentPtrType * _aux;// keep a pointer to maintain assignment of iterator !
  };

 
 // specialisation when no auxiliary data is present 
 template< typename BidirectionalIteratorType, typename Dressing> 
  struct dressed_iterator<BidirectionalIteratorType, Dressing,void> : 
  public boost::iterator_facade<dressed_iterator<BidirectionalIteratorType,Dressing>, Dressing, boost::bidirectional_traversal_tag,Dressing > {
   public : 
    dressed_iterator () {}
    dressed_iterator (BidirectionalIteratorType const & it): _it(it) {}
    BidirectionalIteratorType const & get() const { return _it;}
    BidirectionalIteratorType       & get()       { return _it;}
    operator BidirectionalIteratorType() const { return _it;}
   private:
    friend class boost::iterator_core_access;
    void increment(){ ++_it;}
    void decrement(){ --_it;}
    bool equal(dressed_iterator const & other) const { return(other._it==_it);}
    Dressing dereference() const { return Dressing(_it); } 
    BidirectionalIteratorType _it;
  };


}}
#endif
