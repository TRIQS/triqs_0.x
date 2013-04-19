
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

#ifndef TRACESLICE_H
#define TRACESLICE_H

#include <vector>
#include "hloc.hpp"
#include "small_matrix.hpp"

/**
   One slice in the trace calculation
*/
template<typename VALTYPE>
class TraceSlice { 
  const Hloc & H;
  Hloc::DiagonalOperator & mydiagop; // for temporary calculations
  std::vector<VALTYPE> memChunk;
  std::vector<const Hloc::Bloc *> BlocsOut;
  std::vector<double> Exp_H_tau_acc;
  bool is_nul_;
  std::vector< std::vector < const Hloc::Operator * > >  & NonVanishingOpsOnBlock;
  std::set<const Hloc::Operator *> _nonVanishingOperators;
public:
  // Constructor
  TraceSlice(const Hloc & H_,                       // Hloc
	     Hloc::DiagonalOperator & mydiagop_,    // common workspace, located in the TraceSlice_Stack
	     int mysize,                            // total size of the workspace
	     std::vector< std::vector < const Hloc::Operator * > >  & NonVanishingOpsOnBlock_,// ref : in traceStack
	     bool setAsBoundary = false);
  
  /// Copy Constructor (not expected to be called - no implementation}
  TraceSlice(const TraceSlice &);
  
  /**
     Returns a set of Operator that will not cancel *this
   */
  const std::set<const Hloc::Operator *> & nonVanishingOperators();

  /**
   * Sets the Slice to Op * D * S
   * D is a diagonal operator.
   * If D is NULL, it is interpreted as if D =1 and skipped in the computation 
   */
  TraceSlice & setFrom_Op_D_Slice(const Hloc::Operator &Op, const Hloc::DiagonalOperator *D, const TraceSlice * S);

   /**
   * Sets the Slice to Op * D * S
   * D is a diagonal operator.
   * If D is NULL, it is interpreted as if D =1 and skipped in the computation 
   */
  inline TraceSlice & setFrom_Op_D_Slice1(const Hloc::Operator &Op, const TraceSlice * S);

  /**
   * Sets the Slice to Op * U(delta_tau) * S
   * NB : no check on the sign of dt
   */ 
  TraceSlice & setFrom_Op_U_Slice(const Hloc::Operator &Op, double dt, const TraceSlice * S) {
    assert(S); 
    _nonVanishingOperators.clear();
    // Update the Exp_H_tau_acc from the previous slice (or NULL) if Block is non NULL
    for (int u =0; u<H.NBlocks; ++u)
      if (S->BlocsOut[u] !=NULL)
	Exp_H_tau_acc[u]= S->Exp_H_tau_acc[u] - dt* (S->BlocsOut[u]->H[0] - H.E_GS);
    
    // All blocks are of size 1, no need to compute a diagonal operator, everything is in Exp_H_tau_acc
    if (H.MaxDimBlock ==1)  return setFrom_Op_D_Slice1(Op ,S);
    
    // Computes U(t1,t2) into mydiagop
    for (Hloc::BlocIterator B = H.BlocBegin(); !B.atEnd(); ++B) {
      mydiagop[B][0]=1;
      for (int i =1;  i<B->dim; ++i) mydiagop[B][i] = exp(- dt * B->deltaH[i]); 
    }

    return setFrom_Op_D_Slice(Op,&mydiagop,S);
  }

  /** 
   * Inner product this * D * S
   * D is a diagonal operator.
   * If D is NULL, it is interpreted as if D =1 and skipped in the computation 
   */
  static VALTYPE Slice_D_Slice (const TraceSlice * S1, const Hloc::DiagonalOperator * D, const TraceSlice * S2)  {
    return Slice_D_Slice_internal(S1,D,S2,S1->Exp_H_tau_acc);
  }
  
  /** 
   * Inner product this * U * S
   * D is a diagonal operator.
   * If D is NULL, it is interpreted as if D =1 and skipped in the computation 
   */
  static VALTYPE Slice_U_Slice (const TraceSlice * S1, double dt, const TraceSlice * S2)  {
    std::vector<double> Exp_H_tau_acc_bis(S1->Exp_H_tau_acc.size());

    // Update the Exp_H_tau_acc from the previous slice (or NULL) if Block is non NULL
    for (int u =0; u<S1->H.NBlocks; ++u)
      if (S1->BlocsOut[u] !=NULL)
	Exp_H_tau_acc_bis[u]= S1->Exp_H_tau_acc[u] - dt* (S1->BlocsOut[u]->H[0] - S1->H.E_GS);
    
    // All blocks are of size 1, no need to compute a diagonal operator, everything is in Exp_H_tau_acc
    if (S1->H.MaxDimBlock ==1) 
      return Slice_D_Slice_internal(S1,NULL,S2,Exp_H_tau_acc_bis);
    
    Hloc::DiagonalOperator & mydiagop (const_cast<TraceSlice *>(S1)->mydiagop);
    // Computes U(t1,t2) into mydiagop
    for (Hloc::BlocIterator B = S1->H.BlocBegin(); !B.atEnd(); ++B) {
	mydiagop[B][0]=1;
	for (int i =1;  i<B->dim; ++i) mydiagop[B][i] = exp(- dt * B->deltaH[i]); 
      }

    return Slice_D_Slice_internal(S1, &mydiagop ,S2,Exp_H_tau_acc_bis);
  }

  /// Print
  template<typename TT> 
  friend std::ostream & operator<< (std::ostream & out, const TraceSlice<TT> * S);

  inline bool is_nul() const {return is_nul_;}

private: 

  static VALTYPE Slice_D_Slice_internal (const TraceSlice * S1, const Hloc::DiagonalOperator * D, const TraceSlice * S2,
					 const std::vector<double> & Exp_H_tau_acc_S1);

};


//****************************************************************
// Constructor
template<typename VALTYPE>
TraceSlice<VALTYPE>::TraceSlice(const Hloc & H_, Hloc::DiagonalOperator & mydiagop_,
				int mysize,
				std::vector< std::vector < const Hloc::Operator * > >  & NonVanishingOpsOnBlock_,// ref : in traceStack
				bool setAsBoundary): 
  H(H_), mydiagop(mydiagop_), memChunk(mysize), 
  BlocsOut(H.NBlocks,(Hloc::Bloc*)NULL),
  Exp_H_tau_acc(H.NBlocks,0),is_nul_(false),NonVanishingOpsOnBlock(NonVanishingOpsOnBlock_) {
  if (setAsBoundary) { 
    // I set the matrix as unit
    VALTYPE * restrict pLoc(&memChunk[0]);
    for (Hloc::BlocIterator B = H.BlocBegin(); !B.atEnd(); ++B) {
      BlocsOut[B->num] = B;
      SmallMatrix<VALTYPE,ByLines> M(B->dim,B->dim,pLoc);
      for (int i =0; i<B->dim; ++i) 
	for (int j =0; j<B->dim; ++j) 
	  M(i,j) = (i==j ? 1 : 0);
      pLoc += B->dim*B->dim;
    }
  }
}

//****************************************************************

template<typename VALTYPE>
const std::set<const Hloc::Operator *> & TraceSlice<VALTYPE>::nonVanishingOperators() {

  if ( (_nonVanishingOperators.size() ==0) && (!is_nul())) { //need to compute ite
    for (int u =0; u<H.NBlocks; ++u) {
      if (BlocsOut[u] != NULL) 
	_nonVanishingOperators.insert(NonVanishingOpsOnBlock[ BlocsOut[u]->num ].begin(),
				      NonVanishingOpsOnBlock[ BlocsOut[u]->num ].end());
    }
  }
  return _nonVanishingOperators;
}  

 //****************************************************************

template<typename VALTYPE>
inline TraceSlice<VALTYPE> & TraceSlice<VALTYPE>::setFrom_Op_D_Slice1(const Hloc::Operator & Op, const TraceSlice * S)  {
  assert(S);
  VALTYPE * restrict pLoc(&memChunk[0]);
  const VALTYPE * restrict pS((VALTYPE*)&S->memChunk[0]);
  is_nul_=true;
  for (int u =0; u<H.NBlocks; ++u) {
    const Hloc::Bloc * Bout(S->BlocsOut[u]);
    BlocsOut[u] = NULL;
    if (Bout != NULL) {
      is_nul_=false;
      BlocsOut[u] = Op[Bout].Btarget; 
      if (BlocsOut[u] != NULL) { // Now fill in the Slice
	pLoc[0] = Op[Bout].M.data_()[0] * pS[0];
	pLoc ++;
      }
      pS ++;
    }
  }
  _nonVanishingOperators.clear();
  return *this;
}

 //****************************************************************

template<typename VALTYPE>
TraceSlice<VALTYPE> & TraceSlice<VALTYPE>::setFrom_Op_D_Slice(const Hloc::Operator & Op, 
							      const Hloc::DiagonalOperator * D, const TraceSlice * S)  {
  assert(S);
  VALTYPE * restrict pLoc(&memChunk[0]);
  VALTYPE * restrict pS((VALTYPE*)&S->memChunk[0]);
  is_nul_=true;
  for (Hloc::BlocIterator B = H.BlocBegin(); !B.atEnd(); ++B) {
    const Hloc::Bloc * Bout(S->BlocsOut[B->num]);
    BlocsOut[B->num] = (Bout != NULL ? Op[Bout].Btarget : NULL);
    if (Bout != NULL) {
      is_nul_=false;
      const Hloc::Operator::BlocMatrixElement & OpB(Op[Bout]);
      const int n2 (B->dim), n3 (Bout->dim);
      if (BlocsOut[B->num] != NULL) { // Now fill in the Slice
	int n1 = OpB.M.n1;
	assert(OpB.M.n2 == n3);
	if ( (n1==1) && (n2==1) && (n3==1))
	  pLoc[0] = OpB.M.data_()[0] * pS[0];
	else { 
	  const SmallMatrix<VALTYPE,ByColumns> MSlice(n3,n2,pS);
	  if (D!=NULL)
	    SmallMatrix<VALTYPE,ByColumns>(n1,n2,pLoc).setTo_ADB(OpB.M, (*D)[Bout], MSlice);
	  else
	    SmallMatrix<VALTYPE,ByColumns>(n1,n2,pLoc).setTo_AB(OpB.M, MSlice);
	}
	pLoc += n1*n2;
      }
      pS += n2*n3;
    }
  }
  _nonVanishingOperators.clear();
  return *this;
}

//****************************************************************

/**
 * Do the inner product S1 * U * S2
 * Careful that one assume that the *transpose* of a matrix is in S1
*/
template<typename VALTYPE>
VALTYPE TraceSlice<VALTYPE>::Slice_D_Slice_internal (const TraceSlice * S1,
						     const Hloc::DiagonalOperator * DiagOp,
						     const TraceSlice * S2,
						     const std::vector<double> & Exp_H_tau_acc_S1)  {
  assert(S1); assert(S2); // S1 and S2 have column-ordered matrices
  VALTYPE sum(0);//, sum2(0);
  const VALTYPE * restrict pS1(&S1->memChunk[0]);
  const VALTYPE * restrict pS2(&S2->memChunk[0]);

  for (Hloc::BlocIterator B = S1->H.BlocBegin(); !B.atEnd(); ++B) {
    const Hloc::Bloc * Bout1(S1->BlocsOut[B->num]);
    const Hloc::Bloc * Bout2(S2->BlocsOut[B->num]);
    if (Bout1 == Bout2 && Bout1 != NULL) {   // Bout1 and Bout2 must be identical to contribute
      VALTYPE sum_part = 0;
      if (DiagOp) {  
	for (int i = 0, u=0; i < B->dim; ++i) 
	  for (int k = 0; k < Bout1->dim; ++k, ++u) 
	    sum_part += pS1[u] * ((*DiagOp)[Bout1][k]) * pS2[u];
      }
      else {
	for (int i = 0, u=0; i < B->dim; ++i) 
	  for (int k = 0; k < Bout1->dim; ++k, ++u) 
	      sum_part += pS1[u] * pS2[u];
      }
      sum += sum_part * exp(Exp_H_tau_acc_S1[B->num] + S2->Exp_H_tau_acc[B->num]);
      //sum2 += abs(sum_part * exp(Exp_H_tau_acc_S1[B->num] + S2->Exp_H_tau_acc[B->num]));
     
    }
    if (Bout1 != NULL) pS1 += B->dim*Bout1->dim;
    if (Bout2 != NULL) pS2 += Bout2->dim*B->dim;
  }
  //  assert( abs(abs(sum/sum2) - 1)< 1.e-13); 
  return sum;
}

//*******************************************************************************

template<typename TT> 
std::ostream & operator<< (std::ostream & out, const TraceSlice<TT> * S) {
  if (S == (TraceSlice<TT> *)NULL) return out << "NULL";
    out << "Matrices: " << S->memChunk;
    out << "Connect: " << std::endl;
    for (std::vector<const Hloc::Bloc *>::const_iterator p = S->BlocsOut.begin(); p != S->BlocsOut.end(); ++p) {
      (*p ? out << (*p)->num << " " : out << "NULL ");
    }
    out<<std::endl<<" Exp_H_tau_acc = ";
    std::copy(S->Exp_H_tau_acc.begin(), S->Exp_H_tau_acc.end(), std::ostream_iterator<double>(out,", "));
    out << std::endl;
    return out;
}

#endif
