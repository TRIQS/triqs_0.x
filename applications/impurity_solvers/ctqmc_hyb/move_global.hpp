
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

#ifndef GLOBAL_MOVES_H
#define GLOBAL_MOVES_H

#include "configuration.hpp"
#include "triqs/mc_tools/random_generator.hpp"


class Global_Move { 
  const std::string name;
  Configuration & Config;
  //triqs::mc_tools::random_generator & Random;  
  const std::vector<const Hloc::Operator*> mapping;
  std::vector<Configuration::DET_TYPE*> dets_save;
public :  
  
  typedef std::complex<double> mc_weight_type;

  
  Global_Move(std::string name_, Configuration & Config_, triqs::mc_tools::random_generator & RNG, 
	      const std::vector<const Hloc::Operator*> & mapping_ ):
    name(std::string("Global_Move_") + name_),
    Config(Config_), //Random(RNG),
    mapping(mapping_),
    dets_save(Config.Na,(Configuration::DET_TYPE*)NULL)  {
    // build auxiliary determinants (empty)
    for (int a= 0; a<Config.Na; ++a) {
      Configuration::DET_TYPE * d = Config.dets[a];
      dets_save[a] = new Configuration::DET_TYPE(*(Config_.Delta_tau_proxy[a]), d->size());
    }
  }

  //---------------------

  ~Global_Move() { for (int a= 0; a<Config.Na; ++a) delete dets_save[a];}

  //---------------------

  mc_weight_type attempt() {

#ifdef DEBUG
    std::cout << "I AM IN attempt for Global_Move" << std::endl;
    std::cout << "CONFIG BEFORE: " << Config.DT << std::endl;
#endif
 
    // We first try the change of the DT
    Config.DT.applyGlobalFunction(mapping);
 
    // Now, we will look at the new operators in the trace and reconstruct
    // the other structure of the configuration, like the determinants, Insertions list etc....

    // reserve space to gather the OP_REF on the new operators C, Cdagger.
    std::vector< std::vector<Configuration::OP_REF> > C1(Config.Na),Cdag1(Config.Na),C(Config.Na),Cdag(Config.Na);
    for (int a= 0; a<Config.Na; ++a) {
      C1[a].reserve(100); Cdag1[a].reserve(100);
      C[a].reserve(100); Cdag[a].reserve(100);
    }

    // Go over the list of local operators in DT
    for (Configuration::OP_REF it=Config.DT.OpRef_begin(); ! it.atEnd(); ++it) {
      Configuration::BlockInfo & INFO(Config.info[it->Op->Number]); 
      if (!INFO.isFundamental()) continue;
      if (INFO.dagger) 
	Cdag1[INFO.a].push_back(it);
      else 
	C1[INFO.a].push_back(it);
    }
    
    //-----------------------------------------
    //    Starts recomputation of determinants 
    //-----------------------------------------

    // reverse the order of the C, since they are stored in decreasing time.
    for (int a= 0; a<Config.Na; ++a) {
      for (uint u =0; u<C1[a].size(); ++u) C[a].push_back(C1[a][C1[a].size() - u -1]);
      for (uint u =0; u<Cdag1[a].size(); ++u) Cdag[a].push_back(Cdag1[a][Cdag1[a].size() - u -1]);
    }
    
    // I stored the current det into a copy and 
    // recompute the determinants
    double r1=1, r2=1;
    for (int a= 0; a<Config.Na; ++a) {
      std::swap(dets_save[a],Config.dets[a]);
      Config.dets[a]->recomputeFrom(Cdag[a],C[a]);
      r1 *= Config.dets[a]->determinant();
      r2 *= dets_save[a]->determinant();
    }		
    r1 /=r2;
    //double detratio =  (std::isfinite(real(r1)) ? real(r1) : 0); //+ I* (std::isfinite(imag(r1)) ? imag(r1) : 0);
    // CHANGE Double to the type of determinant to be computed from the class !!!
    double detratio =  (std::isfinite((r1)) ? (r1) : 0); //+ I* (std::isfinite(imag(r1)) ? imag(r1) : 0);

#ifdef DEBUG
    cout<<"About to print"<<endl;
    std::cout << " trace RATIO: " << Config.DT.ratioNewTrace_OldTrace() << endl<<" det ratio"<<detratio << std::endl;
    std::cout << "CONFIG AFTER: " << Config.DT << std::endl;
#endif

    return double(Config.DT.ratioNewTrace_OldTrace()) * detratio;
  }
  
  //-------------------------------------------------------------------------
  
  mc_weight_type accept() { 

    Config.DT.confirm_applyGlobalFunction();

    // clean the map 
    for (Configuration::O_Odag_Insertions_map_type::iterator it = Config.O_Odag_Insertions.begin();
	 it != Config.O_Odag_Insertions.end(); ++it) {
      (*it).second.clear();
    }
    
    // Go over the list of local operators in DT and reconstruct the InsertionTablePtr
    for (Configuration::OP_REF it=Config.DT.OpRef_begin(); ! it.atEnd(); ++it) {
      Configuration::BlockInfo & INFO(Config.info[it->Op->Number]); 
      if (INFO.isFundamental()) continue; 
      //  non fundamental operator
      // if the operator is part of an O_Odag insertion
      if (INFO.InsertionTablePtr) { 
	INFO.InsertionTablePtr->push_back(it);
      }
    }

#ifdef DEBUG
    std::cout << "I AM IN accept for Global_Move" << std::endl;
    std::cout << "CONFIG AFTER: " << Config.DT << std::endl;
#endif
    Config.update_Sign();
    return Config.ratioNewSign_OldSign();
  }   
  
  //-------------------------------------------------------------------------
  
  void reject() {
    
    Config.DT.undo_applyGlobalFunction();
    
#ifdef DEBUG
    std::cout << "I AM IN reject for Global_Move" << std::endl;
    std::cout << "CONFIG AFTER: " << Config.DT << std::endl;
#endif

    for (int a= 0; a<Config.Na; ++a) std::swap(dets_save[a],Config.dets[a]);
  }
  
};

#endif

