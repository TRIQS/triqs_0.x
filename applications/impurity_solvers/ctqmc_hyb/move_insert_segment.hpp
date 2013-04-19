
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

#ifndef MOVE_INSERT_REMOVE_SEGPIC_H_3hurh
#define MOVE_INSERT_REMOVE_SEGPIC_H_3hurh

#include "configuration.hpp"
#include "ctqmc_utils.hpp"
#include <triqs/mc_tools/random_generator.hpp>

/*
  Implementation of the Segment Picture version of the Insert/remove Moves.
  This is only valid if C Cdag alternate for each a level.
 */

/************************

   Insertion of C, C^dagger operator

****************************/
class Insert_Cdag_C_Delta_SegmentPicture { 
  double deltaTau;
  bool no_trivial_reject;
  Configuration & Config;
  triqs::mc_tools::random_generator & Random;
  const int a_level, Nalpha;
  const std::string name;
  triqs::mc_tools::histogram_binned & HISTO_Length_Kinks_Proposed, & HISTO_Length_Kinks_Accepted;
  Configuration::DET_TYPE * det;
  double try_insert_length_max;
public :  
  typedef std::complex<double> mc_weight_type;


  //-----------------------------------------------

  Insert_Cdag_C_Delta_SegmentPicture(int a, Configuration & Config_, triqs::mc_tools::HistogramBinnedMap & HistoMap, triqs::mc_tools::random_generator & RNG  ):
   Config(Config_), Random(RNG),
   a_level(a), Nalpha(Config.COps[a].size()), 
   name( to_string("Insert_Cdagger_C_SegmentPicture_",a)),
   HISTO_Length_Kinks_Proposed (ctqmc_utils::map_insert_check(HistoMap, this->name + "_histo_proposed",triqs::mc_tools::histogram_binned(0,Config.Beta))),
   HISTO_Length_Kinks_Accepted (ctqmc_utils::map_insert_check(HistoMap, this->name + "_histo_accepted",triqs::mc_tools::histogram_binned(0,Config.Beta)))
 {assert(Nalpha==1);}

  //---------------------

  mc_weight_type attempt() {
   det = Config.dets[a_level]; 
   // NB the det pointer has to be recomputed each time, since global moves will change it

#ifdef DEBUG
   std::cout << "I AM IN attempt for Insert_Cdag_C_Delta_SegmentPicture" << std::endl;
   std::cout << "CONFIG BEFORE: " << Config.DT << std::endl;
   // for (int a = 0; a<Config.Na; ++a) print_det(Config.dets[a]);
#endif

   // Pick up a time to insert the first operator
   double tau1 = Random(Config.Beta);

   // Now find the operator A on the same a level at later time, with cyclicity
   // and compute the maximal *oriented* length between the 2 operators (with cyclicity)
   Configuration::OP_REF A;
   if (det->size()>0) { // non empty
    Configuration::DET_TYPE::Cdagger_iterator itCdag = det->Cdagger_begin();
    Configuration::DET_TYPE::C_iterator       itC    = det->C_begin(); 
    while ( (itCdag != det->Cdagger_end()) && (itCdag->tau > tau1) ) {++itCdag;}
    while ( (itC != det->C_end()) &&  (itC->tau > tau1) ) {++itC;}
    if (itCdag != det->Cdagger_begin()) --itCdag; else itCdag = --det->Cdagger_end(); 
    if (itC != det->C_begin()) --itC; else itC = --det->C_end();
    double rC = Config.CyclicOrientedTimeDistance((*itC)->tau - tau1);
    double rCdag = Config.CyclicOrientedTimeDistance((*itCdag)->tau - tau1);
    A = ( rC > rCdag ? *itCdag : * itC);
    try_insert_length_max  = std::min(rC,rCdag);
   }
   else { // empty case.
    try_insert_length_max  = Config.Beta;
    A  = Config.DT.OpRef_end();
   }    

   // pick up the actual segment length
   double rr= Random(try_insert_length_max-2*ctqmc_utils::EPSILON) + ctqmc_utils::EPSILON;

   // deduce the time of the second operator
   double tau2 = Config.CyclicOrientedTimeDistance(tau1 + rr);

   // record the length of the kinks
   if (Config.RecordStatisticConfigurations) {
    deltaTau = Config.CyclicOrientedTimeDistance(tau2 - tau1);
    HISTO_Length_Kinks_Proposed<< abs(deltaTau);
   }

   // Choose the operators to be inserted
   const Hloc::Operator 
    & OpCdag(*Config.CdagOps[a_level][0]),
    & OpC(*Config.COps[a_level][0]);

   // shall we add C^+ C or  C C^+
   // we look whether the next operator is a dagger
   bool Op1_is_dagger = true;
   if (!A.atEnd()) {
    const Configuration::BlockInfo & INFO(Config.info[A->Op->Number]); 
    assert(INFO.isFundamental());
    Op1_is_dagger = INFO.dagger;
   }

   // Insert the operators Op1 and Op2. 
   // O1 will always be the dagger : 
   // Cf doc of insertTwoOperators, order of output OPREF is the same as input operators 
   Configuration::OP_REF O1, O2;
   std::tie (no_trivial_reject,O1,O2) = (Op1_is_dagger ? 
     Config.DT.insertTwoOperators(tau1,OpCdag,tau2,OpC) : 
     Config.DT.insertTwoOperators(tau2,OpCdag,tau1,OpC));
   if (!no_trivial_reject) return 0;

   double tauCdag = (Op1_is_dagger ? tau1 : tau2);
   double tauC = (Op1_is_dagger ? tau2 : tau1);

   // Find the position for insertion in the determinant
   // NB : the determinant store the C in decreasing order.
   int numCdag=1;
   for (Configuration::DET_TYPE::Cdagger_iterator p= det->Cdagger_begin();
     (p != det->Cdagger_end()) && (p->tau > tauCdag)  ; ++p, ++numCdag) {}
   int numC=1;
   for (Configuration::DET_TYPE::C_iterator p= det->C_begin(); 
     (p != det->C_end()) &&  (p->tau > tauC) ; ++p, ++numC) {}

   // acceptance probability
   mc_weight_type p = Config.DT.ratioNewTrace_OldTrace() * det->try_insert(numCdag-1,numC-1,O1,O2);
   int Na(det->size()+1); 
   // !!! det not modified until det->complete_operation is called, so I need to compensate by +1
   double Tratio = Config.Beta * try_insert_length_max / (Na ==1 ? 1 : 2*Na);
   // (Na... term : cf remove move...
#ifdef DEBUG
   std::cout << "Trace Ratio: " << Config.DT.ratioNewTrace_OldTrace() << std::endl;
   std::cout << "p*T: " << p*Tratio << std::endl;
   std::cout << "CONFIG AFTER: " << Config.DT << std::endl;
   //for (int a = 0; a<Config.Na; ++a) print_det(Config.dets[a]);
#endif

   return p*Tratio;
  }

  //----------------

  mc_weight_type accept() { 
   Config.DT.confirm_insertTwoOperators();
   det->complete_operation(); 
   if (Config.RecordStatisticConfigurations) {
    HISTO_Length_Kinks_Accepted << deltaTau;
   }
#ifdef DEBUG
   std::cout << "CONFIG ACCEPT: " << Config.DT << std::endl;
   for (int a = 0; a<Config.Na; ++a) print_det(Config.dets[a]);
#endif
   Config.update_Sign();
   return Config.ratioNewSign_OldSign();
  }

  //----------------

  void reject() {
   if (no_trivial_reject) { 
    Config.DT.undo_insertTwoOperators(); //nothing to be done for the det 
   }
#ifdef DEBUG
   std::cout << "CONFIG REJECT: " << Config.DT << std::endl;
   //for (int a = 0; a<Config.Na; ++a) print_det(Config.dets[a]);
#endif
  }

};

/************************

  Removal of C, C^dagger operator

 ****************************/

class Remove_Cdag_C_Delta_SegmentPicture  { 
 Configuration & Config;
 triqs::mc_tools::random_generator & Random;
 const int a_level, Nalpha;
 Configuration::DET_TYPE * det;
 public :  

 typedef std::complex<double> mc_weight_type;

 const std::string name;

 //----------------------------------

 Remove_Cdag_C_Delta_SegmentPicture(int a, Configuration & Config_, triqs::mc_tools::random_generator & RNG ):
  Config(Config_), Random(RNG),
  a_level(a), Nalpha(Config.COps[a].size()),
  name( to_string("Remove_Cdagger_C_SegmentPicture_",a))
 {}

 //----------------

 mc_weight_type attempt() {

#ifdef DEBUG
  std::cout << "I AM IN attempt for Remove_Cdag_C_Delta_SegmentPicture" << std::endl;
  std::cout << "CONFIG BEFORE: " << Config.DT << std::endl;
  for (int a = 0; a<Config.Na; ++a) print_det(Config.dets[a]);
#endif

  // the det pointer has to be recomputed each time, since global moves will change it
  det = Config.dets[a_level]; 

  // Pick up a couple of C, Cdagger to remove at random
  const int Na(det->size());
  if (Na==0) return 0;

  // I am selecting the n th C or Cdag operator, starting from small time.
  // Because the det orders the tau in decreasing order (cf insert)
  // I need to iterate over the C, Cdag in reverse.
  int n = Random(2*Na); // pick up n
  // the Cs will look like 0 : c1 c2 c1 c2 c1 c2 c1 c2  : beta
  // where c1, c2 are C,Cdag or Cdag,C
  // from n, I now compute the number of the c1,c2 operator (between 0 and N-1)
  // if n = 2p, n1 = p , n2 = p 
  // if n = 2p+1, n1 = p +1 , n2 = p 
  //       e.g. c1[0] c2[0] c1[1] c2[1] c1[2] c2[2] c1[3] c2[3]
  //        n=    0    1     2     3     4     5     6      7
  //eg n=2,p=1              {          }
  //eg n=3,p=1                    {         }
  //eg n=7,p=3        }                                   {    // be careful of the cyclic condition
  int n1 = 1+ (n+1)/2, n2 = 1 + n/2;
  const bool use_cyclicity (n1 > Na);
  if (use_cyclicity) n1 = 1; // cyclic condition

  // decide whether C1 is C or Cdag
  bool Cisfirst = (det->Cdagger_rbegin()->tau > det->C_rbegin()->tau );
  int numC =    (Cisfirst ? n1: n2);
  int numCdag = (Cisfirst ? n2: n1);
#ifdef DEBUG
  cout<<" Cfirst n n1 n2 "<< Cisfirst<< "  "<<n<< "  "<< n1 <<  "  "<< n2 <<endl;
  cout<<numC<< "  "<<numCdag <<endl;
#endif

  // take an iterator on the couple C,Cdag
  Configuration::DET_TYPE::Cdagger_reverse_iterator  itCdag = det->select_reverse_Cdagger( numCdag );
  Configuration::DET_TYPE::C_reverse_iterator           itC = det->select_reverse_C( numC );

  // Remove the operators from the traces
  Config.DT.removeTwoOperators(*itCdag,*itC);

  // first_point is the first of the couple, next_point is the op on a_level *after* the second
  // with cyclicity
  // if cyclicity was not used, the first point is simply the lowest tau
  // otherwise it is the *highest* tau !
  Configuration::OP_REF first_point, next_point;
  const bool C_before_Cdag(((*itC)->tau < (*itCdag)->tau));
  if ((!use_cyclicity && C_before_Cdag) || (use_cyclicity && !C_before_Cdag)) {
   first_point = *itC;
   ++itC;
   next_point = (itC == det->C_rend() ? *det->C_rbegin() : *itC); // cyclicity check
  }
  else {
   first_point = *itCdag;
   ++itCdag;
   next_point = (itCdag == det->Cdagger_rend() ? *det->Cdagger_rbegin() : *itCdag); // cyclicity check
  }

  double length_max =  ( first_point == next_point ? Config.Beta : 
    Config.CyclicOrientedTimeDistance(next_point->tau- first_point->tau) );

  // Acceptance probability
  mc_weight_type p = Config.DT.ratioNewTrace_OldTrace() * det->try_remove(Na-numCdag, Na-numC);
  double Tratio  =  (Na ==1 ? 1 : 2*Na) / (Config.Beta * length_max);      
  // (Na term : because if we have only 1 couple of C Cdagger, n = 0 and 1 will lead to the same couple
  // and this is the only case like this.
#ifdef DEBUG
  cout<< " length_max "<<length_max<<endl;
  std::cout << "RATIO: " << Config.DT.ratioNewTrace_OldTrace() << std::endl;
  std::cout << "CONFIG AFTER: " << Config.DT << std::endl;
#endif
  return p*Tratio;
}

//----------------

mc_weight_type accept() { 
 Config.DT.confirm_removeTwoOperators(); 
 det->complete_operation(); 
#ifdef DEBUG
 std::cout << "CONFIG ACCEPT: " << Config.DT << std::endl;
 for (int a = 0; a<Config.Na; ++a) print_det(Config.dets[a]);
#endif
 Config.update_Sign();
 return Config.ratioNewSign_OldSign();
}   

//----------------

void reject() { 
 Config.DT.undo_removeTwoOperators(); //nothing to be done for the det 
}

};

#endif
