
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

#ifndef TIME_ORDERED_OPERATOR_LIST_H
#define TIME_ORDERED_OPERATOR_LIST_H

#include <map>

/// Type of the time
enum TimeTypeName { MatsubaraContour, KeldyshContour, TripleContour};

/// Explicit type of the time corresponding to TimeTypeName
template<TimeTypeName T> struct TimeType;
template<> struct TimeType<MatsubaraContour> { typedef double type; };
template<> struct TimeType<KeldyshContour>   { typedef std::pair<bool,double> type;};
template<> struct TimeType<TripleContour>    { typedef std::pair<int,double> type;};

/**
   Time_Ordered_Operator_List<timetype> contains a list of time ordered operators
   of type OPERATORTYPE, on one of the 3 possible contours in time.
   The operators are stored with the corresponding time, and some AUXILIARY_INFO (typically a pointer).
   AUXILIARY_INFO must have default construction, copy, and =.

   The class provides :  
      - types : TAUTYPE, the type of the time tau.
      - insertion/removal operators.
      - iteration over operators : e.g. if C is a  Time_Ordered_Operator_List
	for (Time_Ordered_Operator_List::iterator p= C.begin(); p!=C.end(); ++p) {
             p->tau is the time (of type TimeType<timetype>, hence containing possibly the contour position)
             p->Op is a const OPERATORTYPE & to the operator at that time.
	     p->data is the auxiliary stored with the operator. 

 */
template<TimeTypeName mytimetype, typename OPERATORTYPE, typename AUXILIARY_INFO>
class Time_Ordered_Operator_List { 
public:  
  typedef typename TimeType<mytimetype>::type TAUTYPE;
  const TAUTYPE tmin,tmax; 

  /**
   */
  Time_Ordered_Operator_List( TAUTYPE tmin_, TAUTYPE tmax_):tmin(tmin_),tmax(tmax_){}

  /**
     Copy constructor. Makes a deep copy of the data.
  */
  Time_Ordered_Operator_List (const Time_Ordered_Operator_List<mytimetype,OPERATORTYPE,AUXILIARY_INFO> & X):
    tmin(X.tmin),tmax(X.tmax),thelist(X.thelist) {}

protected:
  /// The list of (tau,Op). 
  typedef std::map<TAUTYPE, std::pair<const OPERATORTYPE *,AUXILIARY_INFO> > THELIST_TYPE;
  typedef typename THELIST_TYPE::iterator THELIST_ITERATOR;
  THELIST_TYPE thelist;  
public:
  struct iteratorValue { TAUTYPE tau; const OPERATORTYPE * Op; AUXILIARY_INFO * data; };
  
  class iterator_impl {
    THELIST_ITERATOR p,end;
    iteratorValue VAL;
    inline void convert(){VAL.tau = p->first; 
      VAL.Op = p->second.first;VAL.data = &p->second.second;}
  public:
    iterator_impl(){}
    iterator_impl(THELIST_TYPE & L,bool atend=false){end=L.end();p=(!atend?L.begin():L.end());}
    iterator_impl(THELIST_TYPE & L, THELIST_ITERATOR p_){ end= L.end(); p=p_;} 
    iterator_impl(const iterator_impl& IT):p(IT.p),end(IT.end){}
    iterator_impl& operator=(const iterator_impl& IT2) { p = IT2.p; end = IT2.end; return *this;}
    iterator_impl& operator++() {++p; return *this;}
    iterator_impl& operator--() {--p; return *this;}
    bool operator==(const iterator_impl& IT2) const { return (IT2.p == p);}
    bool operator!=(const iterator_impl& IT2) const { return (IT2.p != p);}
    iteratorValue * operator->(){ convert(); return &(VAL);}
    iteratorValue & operator*() { convert(); return VAL;}
    bool atEnd () const { return (p==end);}
    THELIST_ITERATOR as_simple_map_iterator() const {return p;}
  };

  inline iterator_impl begin() {return iterator_impl(thelist);}
  inline iterator_impl end() {return iterator_impl(thelist,true);}
 
  typedef iterator_impl iterator;// clang is confused otherwise...

  /**
     Insert a new point at (tau,Op), in the time ordered list.
     Returns a pair results : 
       - a bool indicated whether insertion was successfull. 
         Indeed, if 2 operators are exactly at the same tau (at machine precision),
	 insertion will be rejected. [Fix an old pb of version 1]
       - an iterator to the newly inserted operator.
     NB : the AUXILIARY_INFO field associated to the operator is just default-constructed.
   */
  inline std::pair<bool,iterator> insert(TAUTYPE tau, const OPERATORTYPE & Op) { 
    std::pair<THELIST_ITERATOR,bool> res(thelist.insert(std::make_pair(tau,std::make_pair(&Op,AUXILIARY_INFO()))));
    return std::make_pair(res.second, iterator(thelist,res.first));
  }
  
  /**
     Clear the list
  */
  inline void clear() { thelist.clear();}

  /**
     Number of operators
  */
  inline int size() const { return (int)thelist.size();}

  /**
     Removes the operator pointed to by the iterator it.
     Returns : an iterator pointing to the former *right* neighbour of removed operator
               (possibly this->end())
   */
  inline iterator remove(const iterator & it) { 
    iterator it2(it); ++it2; 
    thelist.erase(it.as_simple_map_iterator());
    return it2;
  }

  /**
     Returns an iterator on the first operator Op with Op->tau > tau (strict)
     If cyclic = false, if there isn't such an operator it returns an iterator to end
     If cyclic = true, in this case, it returns begin()
   */
  inline iterator first_operator_after(TAUTYPE tau, bool cyclic=false) { 
    iterator res(thelist,thelist.upper_bound(tau)); // job is done in map already ...
    return (cyclic && res.atEnd() ? this->begin() : res); 
  }
  

};

#endif




