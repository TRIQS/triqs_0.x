
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

#ifndef TRACESLICESTACK_H
#define TRACESLICESTACK_H

#include "trace_slice.hpp"
#include <stack>

/**
   
  The stack for the TraceSlices
  Doc to be written !!!!

*/
template<typename TRACE_SLICE_TYPE >
class TraceSlice_Stack { 

  std::stack<TRACE_SLICE_TYPE *> content;
  const Hloc & H;
  const int slice_size; // size of the slice to be allocated
  Hloc::DiagonalOperator mydiagop; // one DiagonalOperator used by all traceslice
  // For each Block, I store the list of operators that do not vanish on this block
  std::vector< std::vector < const Hloc::Operator * > > NonVanishingOpsOnBlock;

public:

  const TRACE_SLICE_TYPE TraceSliceBoundary; // one special trace slice for the boundaries of the trace

  /// Constructor
  TraceSlice_Stack(const Hloc & H_, int ninit);

  /// Copy Constructor (not expected to be called - no implementation}
  TraceSlice_Stack(const TraceSlice_Stack & SS): 
    H(SS.H),slice_size(SS.slice_size),mydiagop(H),NonVanishingOpsOnBlock(SS.NonVanishingOpsOnBlock),
    TraceSliceBoundary(H,mydiagop,slice_size,NonVanishingOpsOnBlock,true)
  { std::cout << "THEY CALL ME!!!!!!" << std::endl; assert(0); }

  /// Destructor
  ~TraceSlice_Stack() {
    while (!content.empty()) {
      delete content.top();
      content.pop();
    }
  }
  
  /// Pop the uppermost pointer to TraceSlice if any
  inline TRACE_SLICE_TYPE * pop() {
    if (content.size() > 0) {
      TRACE_SLICE_TYPE * res = content.top(); content.pop(); return res;
    }
    else {
      return new TRACE_SLICE_TYPE(H,mydiagop,slice_size,NonVanishingOpsOnBlock);
    }
  }
  
  /// Push a pointer to TraceSlice in the stack
  inline void push(TRACE_SLICE_TYPE * p) { 
    if (p!=NULL) { content.push(p);}
  }
  
};


//****************************************************************


/**
 * Constructor for TraceSlice_Stack
*/
template<typename TRACE_SLICE_TYPE>
TraceSlice_Stack<TRACE_SLICE_TYPE>::TraceSlice_Stack(const Hloc & H_, int ninit) : 
  H(H_),
  slice_size (H.DimHilbertSpace * H.MaxDimBlock),
  mydiagop(H),
  NonVanishingOpsOnBlock(),
  TraceSliceBoundary(H,mydiagop,slice_size,NonVanishingOpsOnBlock,true)
{
  // for each block, build the list of operators that do not cancel it.
  for (Hloc::BlocIterator B= H.BlocBegin(); !B.atEnd(); ++B) { 
    std::vector < const Hloc::Operator * > tmp;
    for (Hloc::OperatorIterator Op = H.OperatorIteratorBegin(); Op != H.OperatorIteratorEnd(); ++Op) { 
      //Op->second is the operator, Cf hloc.hpp
      if (Op->second[B].Btarget !=NULL) tmp.push_back(& Op->second);
    }
    NonVanishingOpsOnBlock.push_back(tmp);
  }
  
  // Push ninit new TraceSlices
  for (int i = 0; i < ninit; ++i) 
    push(new TRACE_SLICE_TYPE(H,mydiagop,slice_size,NonVanishingOpsOnBlock));


}

#endif
