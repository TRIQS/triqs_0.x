
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

#ifndef TRIQS_ARRAYS_RANGE_H
#define TRIQS_ARRAYS_RANGE_H

#include <ostream>

namespace triqs { namespace arrays { 

 class range {
  int first_, last_, step_;
  public:
  enum { from_start = 0, to_end = -1 };

  range( int first__, int last__, int step__=1):first_(first__), last_(last__), step_(step__) {
   /*    ARRAY_PRECHECK((first_ == from_start) || (last_ == to_end) ||
	 (first_ < last_) && (step > 0) ||
	 (first_ > last_) && (step < 0) ||
	 (first_ == last_), (*this) << " is an invalid range.");
	 ARRAY_PRECHECK((last_-first_) % step == 0,
	 (*this) << ": the stride must evenly divide the range");
	 */  }

   range():first_(0),last_(-1),step_(1) {} // i.e. all

  range(const range& r):first_(r.first_), last_(r.last_), step_(r.step_) {}

  int first() const { return first_;}
  int last () const { return last_;}
  int step() const { return step_;}

  range operator+(int shift) const { return range(first_ + shift, last_ + shift, step_); }

  friend inline std::ostream& operator<<(std::ostream& os, const range& range) {
   os << "range(" << range.first() << "," << range.last() << "," << range.step() << ")"; 
   return os;
  }
 };
 
 class ellipsis : public range { 
  public :
   ellipsis( int first__, int last__, int step__=1): range(first__, last__, step__) {}
   ellipsis() : range() {}
 };

}}//namespace triqs::arrays  
#endif 
