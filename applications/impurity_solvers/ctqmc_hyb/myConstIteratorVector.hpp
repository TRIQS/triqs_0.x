
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

#ifndef MYCONSTITERATORVECTOR_H
#define MYCONSTITERATORVECTOR_H

namespace OP_Tools {

  /// 
  template< typename T>
  class myConstIteratorVector {
    typename std::vector<T>::const_iterator p,end;
  public:
    myConstIteratorVector(const std::vector<T> & V ,bool atend=false){ 
      end= V.end(); p=(!atend ? V.begin():end); }
    myConstIteratorVector(const myConstIteratorVector& IT):p(IT.p),end(IT.end){}
    myConstIteratorVector& operator=(const myConstIteratorVector& IT2) { p = IT2.p; end = IT2.end; return *this;}
    myConstIteratorVector& operator++() {++p; return *this;}
    myConstIteratorVector& operator--() {--p; return *this;}
    bool operator==(const myConstIteratorVector& IT2) { return (IT2.p == p);}
    bool operator!=(const myConstIteratorVector& IT2) { return (IT2.p != p);}
    const T * operator->(){ return &(*p);}
    const T & operator*(){ return *p;}
    bool atEnd () const { return (p==end);}
    operator const T  * () {return &(*p);}
  };
}

#endif

