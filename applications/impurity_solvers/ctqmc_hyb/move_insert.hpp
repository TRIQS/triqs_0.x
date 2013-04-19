
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

#ifndef MOVE_INSERT_REMOVE_H_823esfn
#define MOVE_INSERT_REMOVE_H_823esfn

#include "configuration.hpp"
#include "ctqmc_utils.hpp"
#include <triqs/mc_tools/random_generator.hpp>
#include <triqs/mc_tools/histograms.hpp>

/*
  Implementation of the generic version of the Insert/remove Moves.
  This is only valid if C Cdag alternate for each a level.
 */

/************************

   Insertion of C, C^dagger operator

****************************/
class Insert_Cdag_C_Delta { 
  double deltaTau;
  bool no_trivial_reject;
  Configuration & Config;
  triqs::mc_tools::random_generator & Random;
  const int a_level, Nalpha;
  const std::string name;
  triqs::mc_tools::histogram_binned & HISTO_Length_Kinks_Proposed, & HISTO_Length_Kinks_Accepted;
  Configuration::DET_TYPE * det;
public :  
  typedef std::complex<double> mc_weight_type;


  //-----------------------------------------------

  Insert_Cdag_C_Delta(int a,  Configuration & Config_, triqs::mc_tools::HistogramBinnedMap & HistoMap, triqs::mc_tools::random_generator & RNG ):
    Config(Config_), Random(RNG),
    a_level(a), Nalpha(Config.COps[a].size()), 
    name( to_string("Insert_Cdagger_C_",a)),
   HISTO_Length_Kinks_Proposed (ctqmc_utils::map_insert_check(HistoMap, this->name + "_histo_proposed",triqs::mc_tools::histogram_binned(0,Config.Beta))),
   HISTO_Length_Kinks_Accepted (ctqmc_utils::map_insert_check(HistoMap, this->name + "_histo_accepted",triqs::mc_tools::histogram_binned(0,Config.Beta)))
 {}

  //---------------------

  mc_weight_type attempt() {
#ifdef DEBUG
   std::cout << "I AM IN attempt for Insert_Cdag_C_Delta" << std::endl;
   std::cout << "CONFIG BEFORE: " << Config.DT << std::endl;
   // for (int a = 0; a<Config.Na; ++a) print_det(Config.dets[a]);
#endif
   // Pick up the value of alpha and choose the operators
   const Hloc::Operator 
    & Op1(*Config.CdagOps[a_level][Random(Nalpha)]),
    & Op2(*Config.COps[a_level][Random(Nalpha)]);

   // Choice of times.
   double tau1 = Random(Config.Beta), tau2 = Random(Config.Beta);

   // record the length of the kinks
   if (Config.RecordStatisticConfigurations) {
    deltaTau = Config.CyclicOrientedTimeDistance(tau2 - tau1);
    HISTO_Length_Kinks_Proposed<< deltaTau;
   }

   // Insert the operators Op1 and Op2
   Configuration::OP_REF O1, O2;
   std::tie (no_trivial_reject,O1,O2) = Config.DT.insertTwoOperators(tau1,Op1,tau2,Op2);
   if (!no_trivial_reject) return 0;

   // pick up the determinant 
   // NB ; the det has to be recomputed each time, since global moves will change it
   det = Config.dets[a_level]; 
   int numCdag=1;

   // Find the position for insertion in the determinant
   // NB : the determinant store the C in decreasing order.
   for (Configuration::DET_TYPE::Cdagger_iterator p= det->Cdagger_begin();
     (p != det->Cdagger_end()) && (p->tau > tau1)  ; ++p, ++numCdag) {}
   int numC=1;
   for (Configuration::DET_TYPE::C_iterator p= det->C_begin(); 
     (p != det->C_end()) &&  (p->tau > tau2) ; ++p, ++numC) {}

   // acceptance probability
   mc_weight_type p = Config.DT.ratioNewTrace_OldTrace() * det->try_insert(numCdag-1,numC-1,O1,O2);
   double Tratio = std::pow(2*Nalpha* Config.Beta / double(2*(det->size()+1)), 2);

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
   Config.update_Sign();
   return Config.ratioNewSign_OldSign();
  }

  //----------------

  void reject() {
   if (no_trivial_reject) { 
    Config.DT.undo_insertTwoOperators(); //nothing to be done for the det 
   }
  }


};

/************************

  Removal of C, C^dagger operator

 ****************************/

class Remove_Cdag_C_Delta { 
 Configuration & Config;
 triqs::mc_tools::random_generator & Random;
 const int a_level, Nalpha;
 Configuration::DET_TYPE * det;
 public :  

 typedef std::complex<double> mc_weight_type;

 const std::string name;

 //----------------------------------

 Remove_Cdag_C_Delta( int a, Configuration & Config_, triqs::mc_tools::random_generator & RNG  ):
  Config(Config_), Random(RNG),
  a_level(a), Nalpha(Config.COps[a].size()),
  name( to_string("Remove_Cdagger_C_",a))
 {}

 //----------------

 mc_weight_type attempt() {

#ifdef DEBUG
  std::cout << "I AM IN attempt for Remove_Cdag_C_Delta" << std::endl;
  std::cout << "CONFIG BEFORE: " << Config.DT << std::endl;
  //for (int a = 0; a<Config.Na; ++a) print_det(Config.dets[a]);
#endif

  // the det has to be recomputed each time, since global moves will change it
  det = Config.dets[a_level]; 

  // Pick up a couple of C, Cdagger to remove at random
  const int Na(det->size());
  if (Na==0) return 0;
  int numCdag = 1 + Random(Na);
  int numC = 1 + Random(Na);
  // det is in decreasing order.
  numC    = Na - numC  + 1;
  numCdag = Na - numCdag  + 1;

  // Remove the operators from the traces
  // the -Na +1 is to move backward, to compare with V1 ...
  Config.DT.removeTwoOperators(*det->select_Cdagger( numCdag ),
    *det->select_C( numC ));

  // Acceptance probability
  mc_weight_type p = Config.DT.ratioNewTrace_OldTrace() * det->try_remove(numCdag-1,numC-1);
  double Tratio = std::pow(2*Nalpha* Config.Beta / double(2*Na) ,2);

#ifdef DEBUG
  std::cout << "RATIO: " << Config.DT.ratioNewTrace_OldTrace() << std::endl;
  std::cout << "CONFIG AFTER: " << Config.DT << std::endl;
  //for (int a = 0; a<Config.Na; ++a) print_det(Config.dets[a]);
#endif

  return p/Tratio;
 }

 //----------------

 mc_weight_type accept() { 
  Config.DT.confirm_removeTwoOperators(); 
  det->complete_operation(); 
  Config.update_Sign();
  return Config.ratioNewSign_OldSign();
 }   

 //----------------

 void reject() { 
  Config.DT.undo_removeTwoOperators(); //nothing to be done for the det 
 }


};
#endif
