
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

#ifndef DYNAMIC_TRACE_H
#define DYNAMIC_TRACE_H

#include "hloc.hpp"
#include "trace_slice_stack.hpp"
#include "time_ordered_list.hpp"
#include "time_evolution.hpp"

/** 
   Stores a time ordered list of local operators and the trace of their product, with precomputation system for this trace.

  This class :  
    - contains a Time_Ordered_Operator_List to store the (tau,Operator) couples.
    - can insert 1 or 2 operators at some time
    - can remove them
    - precompute the traces or partial traces and can return it.
    - has a undo/confirm mechanism for the last change and 
    return the current trace and the trace before undo.

    The class is templatized on the TimeEvolution, which is a class
    providing the slice * U  operations.
    It is extracted from the class for future developments.

    The class provides 2 importants types : 
    - OP_REF : a Time_Ordered_Operator_List iterator to an operator.
    (cf Time_Ordered_Operator_List documentation).
    - TAUTYPE : the type of times.

*/
template <typename TIME_EVOLUTION>
class DynamicTrace { 
public:
  // type of the matrix elements.
  typedef typename Hloc::REAL_OR_COMPLEX  REAL_OR_COMPLEX;
  // Cf trace_slice.hpp
  typedef TraceSlice<REAL_OR_COMPLEX>  myTraceSlice;

  // structure stored with the operators.
  struct InfoOnNodes { 
    myTraceSlice *L2R_slice,*R2L_slice, *tmp_slice, *tmp2_slice; 
    InfoOnNodes():L2R_slice(NULL),R2L_slice(NULL),tmp_slice(NULL), tmp2_slice(NULL) {} 
  };

  typedef Time_Ordered_Operator_List<MatsubaraContour,Hloc::Operator,InfoOnNodes> Time_Ordered_Operator_List_TYPE;
  typedef typename Time_Ordered_Operator_List_TYPE::iterator OP_REF;
  typedef typename Time_Ordered_Operator_List_TYPE::TAUTYPE TAUTYPE;

  const Hloc & hloc;
protected:
  Time_Ordered_Operator_List_TYPE * OpList, * OpList_save; // list of operators.
  REAL_OR_COMPLEX CurrentTrace, OldTrace; //Current and old value of the trace
  TraceSlice_Stack<myTraceSlice> SliceStack; // Storage of slices
  TIME_EVOLUTION  TimeEvolution; //
  const myTraceSlice * TraceSliceBoundary_ptr; // a special traceslice to handle the boundary
public:

  /**
     Construct a DynamicTrace in Matsubara...
   */
  DynamicTrace(const Hloc & H, double tmax, double tmin=0):
    hloc(H),
    OpList(new Time_Ordered_Operator_List_TYPE(tmin,tmax)), 
    OpList_save(new Time_Ordered_Operator_List_TYPE(tmin,tmax)), 
    SliceStack(H, 100), TimeEvolution(H),
    it1(OpList->begin()), it2(OpList->begin()), it1_bis(it1), it2_bis(it2) {
    lastop = None;
    CurrentTrace = 1;//H.PartitionFunction(tmax-tmin);
    OldTrace = 1;//H.PartitionFunction(tmax-tmin);
    TraceSliceBoundary_ptr = &SliceStack.TraceSliceBoundary;
  }

  /**
     Copy constructor. Makes a deep copy of the data.
   */
  DynamicTrace(const DynamicTrace<TIME_EVOLUTION> & X):
    hloc(X.hloc),
    OpList(new Time_Ordered_Operator_List_TYPE(*X.OpList)), 
    OpList_save(new Time_Ordered_Operator_List_TYPE(*X.OpList_save)), 
    SliceStack(hloc, 100), TimeEvolution(hloc),
    it1(OpList->begin()), it2(OpList->begin()), it1_bis(it1), it2_bis(it2),
    insertedOperators(X.insertedOperators)
  {
    assert(0);// not tested
    assert (X.lastop==None);lastop = None; 
    // lastop = X.lastop;
    // The copy should be made while operations are completed.
    // otherwise it is complicated : I need to identify the it1, it2,e tc...
    // and in practice, copying DynamicTrace is a slow and rare ops which
    // is not done
    CurrentTrace = X.CurrentTrace;
    OldTrace = X.OldTrace;
    TraceSliceBoundary_ptr = &SliceStack.TraceSliceBoundary;
    has_swapped = X.has_swapped;
    // I avoid to copy all the slices, stacks, etc... I simply recompute the slices.
    // it is a bit slower though...
    recomputeTrace_L2R(OpList->end());
    recomputeTrace_R2L(OpList->begin());
  }

private:
  // forbid the = operator
  void operator=(const DynamicTrace<TIME_EVOLUTION> & X){assert(0);}
public:
  
  ///
  ~DynamicTrace(){ delete OpList; delete OpList_save;}

  /// Ratio current value of trace / preceding one
  REAL_OR_COMPLEX ratioNewTrace_OldTrace() const { 
    assert(lastop!=None); return CurrentTrace/OldTrace;}
 
  /// OP_REF to the first element of the list (or to END if it is empty)
  const OP_REF OpRef_begin() const { return OpList->begin();}

  /// OP_REF to the END of the list (as usual in STL, pointing to the next element after the last one)
  const OP_REF OpRef_end() const { return OpList->end();}

  /// Number of operators in the trace
  int Length() const { return OpList->size();}

  /** 
      Given a time tau, returns a set of Operator * which can be inserted
      without canceling the trace *from the left*.
   */
  const std::set<const Hloc::Operator *> & InsertableOperatorAtTime(double tau) {
    // find the trace slice before time tau
    const myTraceSlice * sliL; TAUTYPE tL;
    std::tie(tL,sliL) = L2R_slice_at_left_of (OpList->first_operator_after(tau));
    return const_cast< myTraceSlice *>(sliL)->nonVanishingOperators();
  }

 /* *****************************************************
    
    Progressive insertion of multiple operators.

  *****************************************************/

  /** 
      DOC TO BE WRITTEN 
  */
  std::tuple<bool,OP_REF,OP_REF> insertOperator (TAUTYPE tau, const Hloc::Operator & Op) {
    assert(lastop==None);
    // THIS ROUTINE ONLY MAKE SENSE FOR TAUTYPE is double. PUT A CHECK
    // insert the 2 operators in the list
    // if insert 1 or 2 has a pb, the move will be rejected at the end.
    // no treatement is necessary for Reject 
    
    bool firstcall = (insertedOperators.size()==0);

    // Insert the operator
    bool ok; OP_REF it;
    std::tie(ok,it) = OpList->insert(tau,Op);
    if (!ok)  return std::make_tuple(false,OpRef_end(),OpRef_end()); 
    insertedOperators.push_back(it);

    // Compute the tmp slice of it
    const myTraceSlice * sliL, * sliR; TAUTYPE tR, tL;
    std::tie(tL,sliL) = L2R_slice_at_left_of (it);
    std::tie(tR,sliR) = (firstcall ? R2L_slice_at_right_of(it) : tmp_slice_at_right_of(it) ); // !! this is TMP not, R2L !!
    if (it->data->tmp_slice==NULL) it->data->tmp_slice = SliceStack.pop(); // new slice tmp for it
    TimeEvolution.Op_U_Slice( it->Op, it->tau, tR, sliR, it->data->tmp_slice);

    OP_REF maxit_forward = OpRef_end();
    it1 = it; ++it;

    // try to close the trace
    double trace = TimeEvolution.Slice_U_Slice(sliL, tL, it2->tau, it2->data->tmp_slice);
    if (abs(trace)>1.e-13) { 
      setNewTraceTo(trace);
    }
    else {
    // We recompute the trace from the left starting after it
      std::tie(ok,maxit_forward) = recomputeTrace_R2L_on_tmp_slices(it1->data->tmp_slice,it1->tau,it,OpList->end());
      --maxit_forward;
      assert (!ok); // or the trace would have closed already
    }

    // Now I compute maxit_backward
    OP_REF maxit_backward=OpRef_end();
    if (it1!=OpRef_begin()) { // if Op is the first operator, maxit_backward
      std::tie(ok,maxit_forward) = recomputeTrace_L2R_on_tmp2_slices(it->data->L2R_slice,it->tau,--it1);
      assert (!ok); // or the trace would have closed already
    }
    return std::make_tuple(true,maxit_forward,maxit_backward);
  }
  
  //---------------------------

  /// Undo 
  void undo_insertOperator() { 
    if (lastop==None) return;
    //assert(lastop==Insert2_Step2);
    undo_common();
    for (uint u=0; u<insertedOperators.size(); ++u) 
      remove_and_clean_slices(insertedOperators[u]);  // remove the operator
    insertedOperators.clear();
  }
  
  //---------------------------
  
  /// Confirm 
  void confirm_insertOperator() {
    //assert(lastop==Insert2_Step2);
    confirm_common();
    OP_REF itR = insertedOperators.front();
    OP_REF itL = insertedOperators.back();
    swap_tmp_R2L_slices(itR,itL);
    recomputeTrace_L2R(itL);
    recomputeTrace_R2L(++itL);
    insertedOperators.clear();
  }

  
  /* *****************************************************
    
     Insertion of 1 operator

  *****************************************************/

  /** 
      Inserts 1 operator in the OperatorList and recompute the new trace.
      Returns a pair (ok, ref) where : 
      - ok : true iif the insertion was successfull. cf Time_Ordered_Operator_List
      - ref is an OP_REF to the newly inserted operator
  */
  std::pair<bool,OP_REF> insertOneOperator (TAUTYPE tau1, const Hloc::Operator & OP1) {
    assert(lastop==None);
    // insert the operator in the list
    bool ok; std::tie(ok,it1) = OpList->insert(tau1,OP1);
    if (!ok)  return std::make_pair (false,it1);

    // get the slices around the operator (maybe )
    const myTraceSlice * sliL, * sliR; TAUTYPE tr, tl;
    std::tie(tl,sliL) = L2R_slice_at_left_of  (it1);
    std::tie(tr,sliR) = R2L_slice_at_right_of (it1);
    
    if (it1->data->R2L_slice==NULL) it1->data->R2L_slice = SliceStack.pop(); // new slice
    TimeEvolution.Op_U_Slice( it1->Op, it1->tau, tr,sliR, it1->data->R2L_slice);
    
    setNewTraceTo(TimeEvolution.Slice_U_Slice(sliL, tl, it1->tau, it1->data->R2L_slice));
    lastop = Insert1;
    return std::make_pair(true,it1);
  }

  //---------------------------

  /// Undo the last insertion of 1 operator
  void undo_insertOneOperator() {
    if (lastop==None) return;
    assert(lastop==Insert1);
    undo_common(); 
    remove_and_clean_slices(it1);
  }

  //---------------------------

  /// Confirm the last insertion of 1 operator.
  void confirm_insertOneOperator(){
    assert(lastop==Insert1);
    confirm_common();
    recomputeTrace_L2R(it1);  // recompute 
    ++it1;
    recomputeTrace_R2L(it1);
  }

 /* *****************************************************
     
    Progressive removal of multiple operators.
     
  *****************************************************/
  
protected:
  std::vector <OP_REF> vector_work_removeTwoOperators_In2steps;

public:
  /// First call.
   std::vector<OP_REF> & removeTwoOperators_In2steps_1 (OP_REF OP1, const Hloc::Operator & Op2) {
    assert(lastop==None); 
    it1 = OP1;

    // Compute maxit2
    const myTraceSlice * sliR; TAUTYPE tR;
    std::tie(tR,sliR) = R2L_slice_at_right_of(it1);
    it1_bis = it1; ++it1_bis;  // I keep it1,it2 to the ops to be removed
    bool ok; OP_REF max_it2;
    std::tie(ok,max_it2) = recomputeTrace_R2L_on_tmp_slices(sliR,tR,it1_bis,OpRef_end());
    assert(!ok); // trace MUST vanish with just O1 !

    // Find the operators of type Op2 between OP1 and max_it2
    vector_work_removeTwoOperators_In2steps.clear(); vector_work_removeTwoOperators_In2steps.reserve(10);//check this for perf !?
    for (OP_REF it = it1; it != max_it2; ++it) 
      if (it->Op->Number == Op2.Number) 
	vector_work_removeTwoOperators_In2steps.push_back(it);

    return vector_work_removeTwoOperators_In2steps;
    // a call may not be followed by anything, so I don't set up the lastop
  }

  //---------------------------------------------------------------------

  /// second call : return the delta_tau_max of the reverse move.
  double removeTwoOperators_In2steps_2 (OP_REF OP2) {
    assert(lastop==None); 

    const myTraceSlice * sliL, * sliR; TAUTYPE tR, tL;
    std::tie(tL,sliL) = L2R_slice_at_left_of (it2);
    std::tie(tR,sliR) = tmp_slice_at_right_of(it2); // here is tmp, not R2L
    setNewTraceTo(TimeEvolution.Slice_U_Slice(sliL,tL,tR, sliR));

#ifdef DEBUG
    std::cout << "it1: " << it1->tau << " it2: " << it2->tau << std::endl;
    std::cout << "TL: " << tL << " TR: " << tR << std::endl;
#endif
    
    // Now compute Delta_tau_max for the reverse insert move.
    // first get max_it
    std::tie(tR,sliR) = R2L_slice_at_right_of(it2); // this time it is the R slice !
    bool ok; OP_REF max_it;
    std::tie(ok,max_it) = recomputeTrace_R2L_on_tmp_slices(sliR,tR,OP2,OpRef_end());
    assert(!ok); // trace MUST vanish without O1 !
    
    double Delta_tau_max = (--max_it)->tau - it1->tau;
    assert(Delta_tau_max>=0);
    lastop =Remove2; // the undo and confirm are the same as for RemoveTwoOperators.
    return Delta_tau_max;
  }
  

  /* *****************************************************
    
     Insertion of 2 operators

  *****************************************************/

  /** 
      Inserts 2 operators in the OperatorList and recompute the new trace.
      Returns a 3-tuple (ok, ref1, ref2) where : 
      - ok : true iif the insertion was successfull. cf Time_Ordered_Operator_List
      - ref1, ref2 are OP_REF to the newly inserted operators
  */
  std::tuple<bool,OP_REF,OP_REF> insertTwoOperators (TAUTYPE tau1, const Hloc::Operator & OP1, TAUTYPE tau2, const Hloc::Operator & OP2) {
    assert(lastop==None);
    // insert the 2 operators in the list
    // if insert 1 or 2 has a pb, the move will be rejected at the end.
    // no treatement is necessary for Reject 
    bool ok;
    std::tie(ok,it1) = OpList->insert(tau1,OP1);
    if (!ok)  return std::make_tuple (false,it1,it1);
    std::tie(ok,it2) = OpList->insert(tau2,OP2);
    if (!ok) { OpList->remove(it1); return std::make_tuple (false,it1,it2);}
    
    has_swapped = (it1->tau > it2->tau);
    if (has_swapped) std::swap(it1,it2); // make sure it2 > it1
    
    // recompute the trace from the left between it1 and it2
    const myTraceSlice * sliL, * sliR; TAUTYPE tR, tL;
    std::tie(tL,sliL) = L2R_slice_at_left_of (it2);
    std::tie(tR,sliR) = R2L_slice_at_right_of(it1);
    
    if (recomputeTrace_R2L_on_tmp_slices(sliR,tR,it1,it2).first) 
      setNewTraceTo(TimeEvolution.Slice_U_Slice(sliL,tL,it2->tau,it2->data->tmp_slice));
    else 
      setNewTraceTo(0);

#ifdef DEBUG
    std::cout<< it1->tau<< " " << it2->tau<<" "<<has_swapped<<std::endl;
    //if ((it1->data->tmp_slice)) std::cout<<"SLI1 "<<(myTraceSlice *)(it1->data->tmp_slice)<<std::endl;
    //if ((it2->data->tmp_slice)) std::cout<<"SLIC2 "<<(myTraceSlice *)(it2->data->tmp_slice)<<std::endl;
    //if (sliL) std::cout<<"SliL"<< (myTraceSlice *)(sliL)<<std::endl;
#endif

    lastop = Insert2;
    return (has_swapped ? std::make_tuple (true,it2,it1) : std::make_tuple (true,it1,it2));
  }
   
  //---------------------------
  
  /// Undo the last insertion of 2 operators
  void undo_insertTwoOperators() { 
    if (lastop==None) return;
    assert(lastop==Insert2);
    undo_common();
    remove_and_clean_slices(it1);  // remove the operator
    remove_and_clean_slices(it2);  // remove the operator
  }
  
  //---------------------------
  
  /// Confirm the last insertion of 2 operators.
  void confirm_insertTwoOperators() {
    assert(lastop==Insert2);
    confirm_common();
    swap_tmp_R2L_slices(it1,it2);
    recomputeTrace_L2R(it2);
    recomputeTrace_R2L(++it2);
  }


  /* *****************************************************
     
     Removal of 1 operator
     
  *****************************************************/
  
  /// Removes 1 operator in the OperatorList and recompute the new trace
  void removeOneOperator (OP_REF OP) {
    assert(lastop==None); assert(OP != OpRef_end());
    it1 = OP;
    const myTraceSlice * sliL, * sliR; TAUTYPE tR, tL;
    std::tie(tL,sliL) = L2R_slice_at_left_of (it1);
    std::tie(tR,sliR) = R2L_slice_at_right_of(it1);
    setNewTraceTo(TimeEvolution.Slice_U_Slice(sliL,tL,tR,sliR));
    lastop = Remove1;
  }

  //---------------------------

  /// Undo the last removal of 1 operator
  inline void undo_removeOneOperators(){ 
    if (lastop==None) return;
    assert(lastop==Remove1); undo_common(); }

  //---------------------------

  /// Confirm the last removal of 1 operator
  inline void confirm_removeOneOperators() { 
    assert(lastop==Remove1);
    confirm_common();
    it1 = remove_and_clean_slices(it1); // returns the operator at the left (later time)
    recomputeTrace_R2L(it1);
    if (it1!= OpList->begin()) recomputeTrace_L2R(--it1);
  }
  
  /* *****************************************************
     
     Removal of 2 operators
     
  *****************************************************/
  
  /// Removes 2 operators in the OperatorList and recompute the new trace.
  void removeTwoOperators (OP_REF OP1, OP_REF OP2) {
    assert(lastop==None); assert (OP1!=OP2);
    it1 = OP1; it2 = OP2;
#ifdef DEBUG
    std::cout << "GNA1 " << it1->tau << std::endl;
    std::cout << "GNA2 " << it2->tau << std::endl;
#endif
    if (it1->tau > it2->tau) std::swap(it1,it2); // make sure it2 > it1
    const myTraceSlice * sliL, * sliR; TAUTYPE tR, tL;
    std::tie(tL,sliL) = L2R_slice_at_left_of (it2);
    std::tie(tR,sliR) = R2L_slice_at_right_of(it1);
#ifdef DEBUG
    std::cout << "it1: " << it1->tau << " it2: " << it2->tau << std::endl;
    std::cout << "TL: " << tL << " TR: " << tR << std::endl;
#endif
    it1_bis = it1;  it2_bis = it2; // I keep it1,it2 to the ops to be removed
    ++it1_bis; --it2_bis;
    if (it1_bis->tau > it2_bis->tau) { // there was no operators between it1 and it2 !
      setNewTraceTo(TimeEvolution.Slice_U_Slice(sliL,tL,tR,sliR));
    }
    else {
      if (recomputeTrace_R2L_on_tmp_slices(sliR,tR,it1_bis,it2_bis).first)
	setNewTraceTo(TimeEvolution.Slice_U_Slice(sliL,tL,it2_bis->tau,it2_bis->data->tmp_slice));
      else 
	setNewTraceTo(0);
    }
    lastop =Remove2;
  }
  
  //---------------------------

  /// Undo the last removal of 2 operators
  inline void undo_removeTwoOperators() {
    if (lastop==None) return;
    assert(lastop==Remove2);undo_common(); 
  }

  //---------------------------

  /// Confirm the last removal of 2 operators.
  void confirm_removeTwoOperators() {
    assert(lastop==Remove2);
    confirm_common();
    if (!(it1_bis->tau > it2_bis->tau)) {
      swap_tmp_R2L_slices(it1_bis,it2_bis);
    }
    // special case where the removed operators are next to each other
    // in this case it1_bis == it2 > it2_bis == it1
    else {
      // there is an operator to the right of both operators
      if (it2_bis != OpList->begin()) {
        --it2_bis;
      }
      // worst situation: both operators are next to the right border
      else {
        ++it1_bis;
        remove_and_clean_slices(it2); 
        remove_and_clean_slices(it1);
        recomputeTrace_R2L(it1_bis);
        return;
      }
    }
    remove_and_clean_slices(it2); 
    remove_and_clean_slices(it1);
    recomputeTrace_L2R(it2_bis);
    recomputeTrace_R2L(++it2_bis);
  }
  
 
 /* *****************************************************
     
     Remove one operator and add another one at a different time
     
  *****************************************************/
 
  /// Removes 1 operator and insert a new one in the OperatorList and recompute the new trace.
  std::pair<bool,OP_REF> insert_and_remove_One_Operator (OP_REF OP_to_remove, TAUTYPE tau, const Hloc::Operator & OP_to_insert) {
    assert(lastop==None);
    it1 = OP_to_remove; 
    bool ok;
    std::tie(ok,it2) = OpList->insert(tau,OP_to_insert);
    if (!ok)  return std::make_pair (false,it1);

    has_swapped = (it1->tau > it2->tau);
    if (has_swapped) std::swap(it1,it2); // make sure it1 <= it2

    const myTraceSlice * sliL, * sliR; TAUTYPE tR, tL;
    std::tie(tL,sliL) = L2R_slice_at_left_of (it2);
    std::tie(tR,sliR) = R2L_slice_at_right_of(it1);
    
    it1_bis = it1;  it2_bis = it2;
    if (has_swapped) --it2_bis; else ++it1_bis; // no possible crossing
    assert(!(it2_bis->tau < it1_bis->tau));

    if (recomputeTrace_R2L_on_tmp_slices(sliR,tR,it1_bis,it2_bis).first)
      setNewTraceTo(TimeEvolution.Slice_U_Slice(sliL,tL,it2_bis->tau,it2_bis->data->tmp_slice));
    else
      setNewTraceTo(0);

    lastop =Insert_Remove1;
    return std::make_pair (true,(has_swapped ? it1 : it2));
  }
  
  //---------------------------

  /// Undo 
  inline void undo_insert_and_remove_One_Operator () {
    if (lastop==None) return;
    assert(lastop==Insert_Remove1);
    undo_common();
    remove_and_clean_slices(has_swapped ? it1 : it2);
  }

  //---------------------------

  /// Confirm 
  void confirm_insert_and_remove_One_Operator () {
    assert(lastop==Insert_Remove1);
    confirm_common();
    remove_and_clean_slices(has_swapped ? it2 : it1);
    swap_tmp_R2L_slices(it1_bis,it2_bis);
    recomputeTrace_L2R(it2_bis);
    recomputeTrace_R2L(++it2_bis);
  }
 

  /* *****************************************************
     
     ApplyGlobalFunction
     
  *****************************************************/
  
  /** 
      Removes 2 operators in the OperatorList and recompute the new trace.
      F is a function : operator_number -> operator 
   */
  void applyGlobalFunction (std::vector<const Hloc::Operator*> const & F){
    lastop = ApplyGlobalFunction;
    if (OpList->size()==0) {OldTrace = CurrentTrace; return;}
    assert (OpList_save->size()==0); // cleaned by accept and reject
    
    for (OP_REF p = OpList->begin(); !p.atEnd() ; ++p) { 
      /*
	uint n = p->Op->Number;
	assert ( n < F.size() );
	const Hloc::Operator * OP = F[p->Op->Number];
	bool ok = OpList_save->insert(p->tau, * OP).first; 
      */
#ifdef NDEBUG
     OpList_save->insert(p->tau, * (F[p->Op->Number])); 
#else
     bool ok = OpList_save->insert(p->tau, * (F[p->Op->Number])).first; 
     assert (ok); 
#endif
    }
    
    std::swap(OpList,OpList_save); // save the operator list

    recomputeTrace_R2L(OpList->begin());
    OP_REF it(OpList->end()); --it; // list is not empty
    setNewTraceTo(TimeEvolution.Slice_U_Slice(TraceSliceBoundary_ptr, OpList->tmax, it->tau, it->data->R2L_slice));
  }
  
  //---------------------------

  /// Undo the last global move
  inline void undo_applyGlobalFunction() { 
    if (lastop==None) return;
    undo_common();
    if (OpList->size()==0) return;
    std::swap(OpList,OpList_save); // restore the original operator list
    // clean the OpList_save
    for (OP_REF p = OpList_save->begin(); p!= OpList_save->end(); ++p) clean_slices(p);
    OpList_save->clear();
    assert (OpList_save->size()==0); // cleaned by accept and reject
  }
  
  //---------------------------
  
  /// Confirm the last global move.
  void confirm_applyGlobalFunction() {
    confirm_common();
    if (OpList->size()==0) return;
    //    assert (OpList->size()!=0);
    // clean the OpList_save
    for (OP_REF p = OpList_save->begin(); p!= OpList_save->end(); ++p) clean_slices(p);
    OpList_save->clear();
    OP_REF it(OpList->end()); --it; // list is not empty
    recomputeTrace_L2R(it);
    assert (OpList_save->size()==0); // cleaned by accept and reject
  }
  

 /* *****************************************************

    Compute the trace with the operator OP replaced by Replacement

  *****************************************************/

  REAL_OR_COMPLEX traceRatioWithOneOperatorReplaced (OP_REF OP, const Hloc::Operator & Replacement) const {
    assert(lastop==None);
    // I overrule the const for TSS. The operation at the end will NOT change TSS
    TraceSlice_Stack<myTraceSlice>  * TSS = const_cast<TraceSlice_Stack<myTraceSlice> *> (&SliceStack);
    // new slice for temporary calculation
    myTraceSlice * tmp = TSS->pop(); 
    // get the slices around the operator
    const myTraceSlice * sliL, * sliR; TAUTYPE tr, tl;
    std::tie(tl,sliL) = L2R_slice_at_left_of  (OP);
    std::tie(tr,sliR) = R2L_slice_at_right_of (OP);
    // Compute sliL Replacement sliR, via tmp
    TimeEvolution.Op_U_Slice( &Replacement, OP->tau, tr,sliR, tmp);
    REAL_OR_COMPLEX res = TimeEvolution.Slice_U_Slice(sliL, tl, OP->tau, tmp);
    TSS->push(tmp); // give back the temporary slice
    assert(std::isfinite(CurrentTrace));
    return res/CurrentTrace;
  }

  //---------------------------------------------------

  ///
  friend std::ostream & operator<< (std::ostream & out, const DynamicTrace & DT) {
    out<<"-------------------------"<<std::endl;
    out<<"Time Ordered List of Operators"<<std::endl;
    for (OP_REF p= DT.OpList->begin(); p!=DT.OpList->end(); ++p) 
      out<<" time = "<<p->tau<<"  Operator "<<p->Op->name<<std::endl; // add aux data if needed...
    out << "Current trace: " << REAL_OR_COMPLEX(DT.CurrentTrace) << std::endl;
    out << "Old trace: " << REAL_OR_COMPLEX(DT.OldTrace) << std::endl;
    out<<"-------------------------"<<std::endl;
    return out;
  }

  inline std::ostream & operator<<(std::ostream & out) {
    out<<"Time Ordered List of Operator"<<std::endl;
    for (OP_REF p= OpList->begin(); p!=OpList->end(); ++p) 
      out<<" time = "<<p->tau<<"  Operator "<<p->Op->name<<std::endl; // add aux data if needed...
    out<<"-------------------------"<<std::endl;
    return out;
  }
 
protected:

  enum LastOperation {None, Insert1, Insert2, Remove1, Remove2, ApplyGlobalFunction,Insert_Remove1};
  LastOperation lastop;
  bool has_swapped;
  OP_REF it1,it2, it1_bis, it2_bis;
  std::vector<OP_REF> insertedOperators;

  inline void setNewTraceTo(REAL_OR_COMPLEX NT) { 
    OldTrace = CurrentTrace;CurrentTrace = NT;
  }
  //-------------------------------------------------

  // TO BE USED IN ALL UNDO/CONFIRM
  inline void undo_common() {CurrentTrace = OldTrace;lastop = None;}
  inline void confirm_common() { lastop = None; }
  
  //---------------------------------------

  // Clean slices around an operator
  inline void clean_slices(OP_REF it) { 
    SliceStack.push(it->data->L2R_slice); // no pb if it is NULL. It does nothing, Cf SliceStack
    SliceStack.push(it->data->R2L_slice);
    SliceStack.push(it->data->tmp_slice);
  }
  //---------------------------------------

  // remove the operator while cleaning the slices used...
  inline OP_REF remove_and_clean_slices(OP_REF it) { 
    clean_slices(it);
    return OpList->remove(it);
  }

  //---------------------------------------

  // returns the R_slice of the operator at the right of it or TraceSliceBoundary_ptr if it is at the boundary
  inline std::pair<TAUTYPE, const myTraceSlice *> R2L_slice_at_right_of (OP_REF it) const { 
    if (it==OpList->begin()) return std::make_pair(OpList->tmin,TraceSliceBoundary_ptr);
    --it; return std::make_pair(it->tau,it->data->R2L_slice); 
  }
  
  //---------------------------------------
  
  // returns the tmp_slice of the operator at the right of it or TraceSliceBoundary_ptr if it is at the boundary
  inline std::pair<TAUTYPE, const myTraceSlice *> tmp_slice_at_right_of (OP_REF it) const { 
    if (it==OpList->begin()) return std::make_pair(OpList->tmin,TraceSliceBoundary_ptr);
    --it; return std::make_pair(it->tau,it->data->tmp_slice); 
  }
  
  //---------------------------------------
  
  // returns the L_slice of the operator at the left of it or TraceSliceBoundary_ptr if it is at the boundary
  inline std::pair<TAUTYPE, const myTraceSlice *> L2R_slice_at_left_of (OP_REF it) const { 
    ++it;
    return (it==OpList->end() ? std::make_pair(OpList->tmax,TraceSliceBoundary_ptr) : 
	    std::make_pair(it->tau, (const myTraceSlice *)it->data->L2R_slice));
  }
  
  //---------------------------------------

  inline void ComputeSlice(OP_REF it, TAUTYPE tstart, const myTraceSlice * sli_start, myTraceSlice * sli)  {
    if (sli==NULL) sli = SliceStack.pop(); 
    TimeEvolution.Op_U_Slice( it->Op, it->tau, tstart, sli_start, sli);
  }

  //---------------------------------------

  // recompute R2L_slice from it ->left boundary including it
  void recomputeTrace_R2L(OP_REF it) { 
    if (it == OpList->end()) return;
    const myTraceSlice * sli; TAUTYPE t_r;
    std::tie(t_r,sli) = R2L_slice_at_right_of(it);
    for (; it != OpList->end(); ++it) {
      if (it->data->R2L_slice==NULL) it->data->R2L_slice = SliceStack.pop();
      TimeEvolution.Op_U_Slice(it->Op, it->tau, t_r, sli ,it->data->R2L_slice);
      //ComputeSlice(it,t_r,sli,it->data->R2L_slice);
      t_r = it->tau;
      sli= it->data->R2L_slice;
    }
  }

  //---------------------------------------

  // recompute L_slice from iter1 to iter2 *included*, starting from the slice slice_start at time t_start
  // while saving the slices for future undo
  // return true is the last slice is NOT null, false otherwise
  // the computation STOPS before iter2 if the trace is all 0
  std::pair<bool,OP_REF> recomputeTrace_R2L_on_tmp_slices(const myTraceSlice * sli_start, TAUTYPE t_start, OP_REF iter1, OP_REF iter2) { 
    //assert (iter2 !=OpList->end()); assert (iter2->tau >= iter1->tau); ++iter2;
    if (iter2 !=OpList->end()) {  assert (iter2->tau >= iter1->tau); ++iter2;}
    OP_REF it;
    for (it = iter1; (it != iter2) && (!sli_start->is_nul()); ++it) {
      if (it->data->tmp_slice==NULL) it->data->tmp_slice = SliceStack.pop();
      TimeEvolution.Op_U_Slice(it->Op, it->tau, t_start ,sli_start, it->data->tmp_slice);
      t_start = it->tau;
      sli_start = it->data->tmp_slice;
    }
    return std::make_pair(!sli_start->is_nul(),it);
  }
 

  //---------------------------------------

  // recompute L2R_slice it -> right boundary it including it
  void recomputeTrace_L2R( OP_REF it) { 
    if (OpList->size()==0) return ;
    const myTraceSlice * sli; TAUTYPE t_l;
    if (it==OpList->end()) { --it;}
    bool fini= false;
    std::tie(t_l,sli) = L2R_slice_at_left_of(it);
    do {
      if (it->data->L2R_slice==NULL) it->data->L2R_slice = SliceStack.pop();
      TimeEvolution.Op_U_Slice( &it->Op->Transpose(), t_l, it->tau, sli, it->data->L2R_slice);
      sli = it->data->L2R_slice;
      t_l = it->tau;
      fini = (it ==OpList->begin());
      if (!fini) --it;
    }
    while  (!fini); 
  }
  
  //---------------------------------------

  // recompute L2R_slice it -> right boundary it including it
  std::pair<bool,OP_REF> recomputeTrace_L2R_on_tmp2_slices(const myTraceSlice * sli_start, TAUTYPE t_start, OP_REF it) { 
    if (OpList->size()==0) return std::make_pair(false,it);
    if (it==OpList->end()) { --it;}
    bool fini= false;
    do {
      if (it->data->tmp2_slice==NULL) it->data->tmp2_slice = SliceStack.pop();
      TimeEvolution.Op_U_Slice( &it->Op->Transpose(), t_start, it->tau, sli_start, it->data->L2R_slice);
      sli_start = it->data->tmp2_slice;
      t_start = it->tau;
      fini = (it ==OpList->begin()) ||  (sli_start->is_nul());
      if (!fini) --it;
    }
    while  (!fini); 
    return std::make_pair(!sli_start->is_nul(),it);
  }

  //---------------------------------------

  // swaps the R and tmp slice for all operators between iter1 and iter2, including them.
  void swap_tmp_R2L_slices(OP_REF iter1, OP_REF iter2) { 
    assert (iter2 !=OpList->end()); 
    assert (iter1 !=OpList->end()); 
    assert (!(iter1->tau > iter2->tau));
    ++iter2;
    for (OP_REF it = iter1; it != iter2; ++it) {
      assert(it->data->tmp_slice);
      std::swap(it->data->tmp_slice,it->data->R2L_slice); 
    }
  }

}; 
#endif
