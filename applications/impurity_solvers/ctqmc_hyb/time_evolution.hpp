
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

#ifndef TIME_EVOLUTION_H
#define TIME_EVOLUTION_H
#include "hloc.hpp"

/**

  The explicit time evolution operator U

*/
template<typename REAL_OR_COMPLEX> class TimeEvolutionSimpleExp {  
  const Hloc & H;

public:
  typedef TraceSlice<REAL_OR_COMPLEX> myTraceSlice;

  TimeEvolutionSimpleExp(const Hloc &h): H(h) {}
   
  inline void Op_U_Slice(const Hloc::Operator * Op, double t1, double t2, 
			 const myTraceSlice * slice_in, myTraceSlice * slice_res ) const { 
    assert (t1>=t2);
    slice_res->setFrom_Op_U_Slice(*Op, t1-t2, slice_in);
  }

  /// returns the inner product slice1 * U(t1,t2) * slice2
  inline Hloc::REAL_OR_COMPLEX Slice_U_Slice( const myTraceSlice * slice1,  double t1, double t2, 
					      const myTraceSlice * slice2) const { 
    assert (t1>=t2);
    return myTraceSlice::Slice_U_Slice (slice1,t1-t2, slice2);
  }

};

#endif
