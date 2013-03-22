
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

#ifndef MOVE_MOVE_CC_H_w3hbv45
#define MOVE_MOVE_CC_H_w3hbv45

#include "configuration.hpp"
#include "triqs/mc_tools/random_generator.hpp"

/************************

   Move one operator

****************************/
class Move_C_Delta { 
  const std::string name;
  Configuration & Config;
  triqs::mc_tools::random_generator & Random;
  Configuration::DET_TYPE * det;
  Configuration::DET_TYPE::RollDirection roll_matrix;

  //-----------------------------------------------

public:

  typedef std::complex<double> mc_weight_type;


  Move_C_Delta(Configuration & Config_, triqs::mc_tools::random_generator & RNG ): name("Move_C_Delta"), 
    Config(Config_), Random(RNG) {}

  //----------------
  
  mc_weight_type attempt() {
    
#ifdef DEBUG
    std::cout << "I AM IN attempt for Move" << std::endl;
    std::cout << "CONFIG BEFORE: " << Config.DT << std::endl;
    for (int a = 0; a < Config.Na; ++a) print_det(Config.dets[a]);
#endif

    // find an operator (anyone) to be moved.
    // to do this, I go over the DT trace
    const int L(Config.DT.Length());
    if (L==0) return 0;
    const int N(Random(L));

    // I go over the trace to find the op oldOpref
    // and then I get its info
    Configuration::OP_REF oldOpref = Config.DT.OpRef_begin();
    for (int i=0; i<N; ++i, ++oldOpref) {}
    int a = Config.info[oldOpref->Op->Number].a;
    int Nalpha = Config.COps[a].size();
    bool isdagger = Config.info[oldOpref->Op->Number].dagger;
    double oldtau = oldOpref->tau;

    // chose the new alpha (this is done here for compatibility)
    int newalpha = Random(Nalpha);

    int num;
    double tR, tL;
    det = Config.dets[a];
    Configuration::DET_TYPE::Cdagger_iterator itCdag = det->Cdagger_begin();
    Configuration::DET_TYPE::C_iterator       itC    = det->C_begin();

    if (det->size() > 1) {

      // find the C and Cdag operators at the right of oldOpref (at smaller times)
      // it might be Cdagger_end() or C_end()
      // also get the number num of oldOpref in the det
      num = 0;
      while ( (itCdag != det->Cdagger_end()) && (itCdag->tau >= oldtau) ) {++itCdag; if(isdagger) ++num;}
      while ( (itC != det->C_end()) &&  (itC->tau >= oldtau) ) {++itC; if(!isdagger) ++num;}

      // find the times of the operator at the right of oldOpref with cyclicity
      // then deduce the closest one and put its distance to oldOpref in tR
      double tRdag = (itCdag != det->Cdagger_end() ? (*itCdag)->tau : det->Cdagger_begin()->tau);
      double tRnodag = (itC != det->C_end() ? (*itC)->tau : det->C_begin()->tau);
      tR = std::min(Config.CyclicOrientedTimeDistance(oldtau - tRdag), Config.CyclicOrientedTimeDistance(oldtau - tRnodag));

      // move itCdag and itC to the operators on the left with cyclicity
      // find their times and deduce the closest one with distance tL
      if (isdagger) --itCdag; else --itC;
      if (itCdag != det->Cdagger_begin()) --itCdag; else itCdag = --det->Cdagger_end();
      if (itC != det->C_begin()) --itC; else itC = --det->C_end();
      tL = std::min(Config.CyclicOrientedTimeDistance((*itCdag)->tau - oldtau), Config.CyclicOrientedTimeDistance((*itC)->tau - oldtau));

    } else {

      num = 1;
      tR = (isdagger ?  Config.CyclicOrientedTimeDistance(oldtau - (*itC)->tau) : Config.CyclicOrientedTimeDistance(oldtau - (*itCdag)->tau));
      tL = (isdagger ?  Config.CyclicOrientedTimeDistance((*itC)->tau - oldtau) : Config.CyclicOrientedTimeDistance((*itCdag)->tau - oldtau));

    }

    // the new operator newOp with a new chosen alpha
    double newtau =  Config.CyclicOrientedTimeDistance(oldtau - tR + Random(tL + tR - 2*ctqmc_utils::EPSILON) + ctqmc_utils::EPSILON);
    const Hloc::Operator * newOp;
    newOp = (isdagger ? Config.CdagOps[a][newalpha] : Config.COps[a][newalpha]);

    // remove oldOpref and insert newOp in the Trace
    bool ok;
    Configuration::OP_REF Op;
    std::tie(ok,Op) = Config.DT.insert_and_remove_One_Operator(oldOpref, newtau, *newOp);
    if (!ok) return 0;

    // in the following we want to see if we need to roll the det
    roll_matrix = Configuration::DET_TYPE::None;

    // check if we went through \tau = 0 or \tau = \beta
    if (newtau - oldtau > tL) roll_matrix = (isdagger ? Configuration::DET_TYPE::Down  : Configuration::DET_TYPE::Right);
    if (oldtau - newtau > tR) roll_matrix = (isdagger ? Configuration::DET_TYPE::Up : Configuration::DET_TYPE::Left);

    // acceptance probability
    mc_weight_type p = Config.DT.ratioNewTrace_OldTrace() * (isdagger ? det->try_change_row(num-1,Op) : det->try_change_col(num-1,Op));

#ifdef DEBUG
    std::cout << "oldtau, tR, tL, newtau: " << oldtau << " " << tR << " " << tL << " " << newtau << endl;
    std::cout << "Trace Ratio: " << Config.DT.ratioNewTrace_OldTrace() << std::endl;
    std::cout << "Det Ratio: " << p/mc_weight_type(Config.DT.ratioNewTrace_OldTrace()) << std::endl;
    std::cout << "Rate: " << p << std::endl;
    std::cout << "CONFIG AFTER: " << Config.DT << std::endl;
    for (int a = 0; a < Config.Na; ++a) print_det(Config.dets[a]);
#endif

    return p;
  }

  //----------------
  
  mc_weight_type accept() { 
    Config.DT.confirm_insert_and_remove_One_Operator();
    det->complete_operation(); 

#ifdef DEBUG
    std::cout << "CONFIG ACCEPT: " << Config.DT << std::endl;
    std::cout << "ROLL MATRIX: " << roll_matrix << std::endl;
    for (int a = 0; a < Config.Na; ++a) print_det(Config.dets[a]);
#endif

    Config.update_Sign();
    return Config.ratioNewSign_OldSign() * det->roll_matrix(roll_matrix);
  }   
  
  //----------------
  
  void reject() { 
    Config.DT.undo_insert_and_remove_One_Operator();

#ifdef DEBUG
    std::cout << "CONFIG REJECT: " << Config.DT << std::endl;
    for (int a = 0; a < Config.Na; ++a) print_det(Config.dets[a]);
#endif

  }
  
};

#endif
