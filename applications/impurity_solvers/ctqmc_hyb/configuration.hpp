
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
#ifndef CONFIGURATION_H_37j3hrf
#define CONFIGURATION_H_37j3hrf

#include "improved_python_dict.hpp"
#include <triqs/gf/block.hpp>
#include <triqs/gf/imtime.hpp>
#include "dynamic_trace.hpp"
#include <triqs/mc_tools/mc_generic.hpp>
#include <map>
#include "hloc.hpp"
#include "gf_binner.hpp"
#include "det_manip_hyb.hpp"

/**
  The configuration of the Monte Carlo
  */
struct Configuration { 
 typedef DynamicTrace< TimeEvolutionSimpleExp <Hloc::REAL_OR_COMPLEX> > DYNAMIC_TRACE;
 typedef DYNAMIC_TRACE::OP_REF OP_REF;

 /**
   Stores all additionnal informations on an operator :
   - for a fundamental operator (a,alpha,dagger)
   - for composite operators, a = -1, and possibly the pointeur to its insertion table
   */
 struct BlockInfo { 
  int a, alpha; bool dagger; 
  std::vector<OP_REF> * InsertionTablePtr; // for the O_Odag_insertion
  BlockInfo():a(-1), alpha(-1), dagger(false),InsertionTablePtr(NULL){} 
  BlockInfo(const BlockInfo & BI):a(BI.a), alpha(BI.alpha), dagger(BI.dagger),InsertionTablePtr(BI.InsertionTablePtr){} 
  bool isFundamental() const { return (a>=0); }
 };

 /// Proxy to call Delta with 2 OP_REF
 struct Delta_Proxy { 
  typedef double result_type;
  typedef OP_REF argument_type;
  const triqs::gf::gf_view<triqs::gf::imtime> & Delta;
  gf_grid_evaluator<triqs::gf::gf_view<triqs::gf::imtime>> Delta_eval; 
  const std::vector<BlockInfo> & info;
  Delta_Proxy (const triqs::gf::gf_view<triqs::gf::imtime> & Delta_, const std::vector<BlockInfo> & info_ ):
   Delta(Delta_),Delta_eval(Delta_), 
   info(info_) { }
  Delta_Proxy (Delta_Proxy const & ) = default;
  //Delta_Proxy (Delta_Proxy && ) = default;
  Delta_Proxy & operator =(Delta_Proxy const & ) = delete;
  double operator()(OP_REF A, OP_REF B) const { 
   return Delta_eval(info[A->Op->Number].alpha,A->tau, 
     info[B->Op->Number].alpha,B->tau);}
 };

 typedef detManip<Delta_Proxy> DET_TYPE;

 /// Object for storing the multiple expansion 
 struct O_Odag_Insertions_type { 
  std::vector<OP_REF> Insertions1, Insertions2;
  void clear(){Insertions1.clear();Insertions2.clear(); Insertions1.reserve(100);  Insertions2.reserve(100);}
 };
 // Stores the relation between a pair of operator and Insertions tables. MUST be initialized by bext routine.
 typedef std::map<std::pair<std::string,std::string>, O_Odag_Insertions_type > O_Odag_Insertions_map_type; 
 O_Odag_Insertions_map_type O_Odag_Insertions; 
 // Creates an entry in this map and updates the InsertionTablePtr in info.
 O_Odag_Insertions_type & create_O_Odag_Insertion (const Hloc::Operator & Op1, const Hloc::Operator & Op2);

 /// Storing the O1(tau) K(tau,tau') O2(tau') insertion
 struct O1_K_O2_Insertions_type { 
  std::vector<std::pair<OP_REF,OP_REF> > Insertions;
 };
 // map will be initialized by the constructor of the moves.
 std::map<std::string, O_Odag_Insertions_type > O1_K_O2_Insertions; 

 ///
 Configuration(triqs::python_tools::improved_python_dict params, Hloc * hloc, triqs::gf::gf_view<triqs::gf::block_index,triqs::gf::gf<triqs::gf::imtime>> &);

 /// 
 // Configuration (const Configuration & C);

 ///
 ~Configuration();

 /// Distance between 2 points ....
 inline double DeltaTau(double t1, double t2) { double r = t2 - t1; return (r>=0 ? r : r+ Beta);}

 /// puts back the time in [0,Beta]
 inline double CyclicOrientedTimeDistance(double x){ double r= fmod(x,Beta); return (r >=0 ? r : r+ Beta);}

 triqs::gf::gf_view<triqs::gf::block_index, triqs::gf::gf<triqs::gf::imtime>> Delta_tau;
 const Hloc & H;
 const double Beta;
 DYNAMIC_TRACE DT;
 const int Na;
 std::vector< DET_TYPE *> dets;  
 std::vector< std::vector<const Hloc::Operator *> > COps, CdagOps; // the operators to be inserted grouped in "a" blocks
 std::vector<Delta_Proxy *> Delta_tau_proxy;
 std::vector<BlockInfo> info; // info[number_of_an_op] gives access to its a, alpha, dagger 
 const bool RecordStatisticConfigurations;

 int ratioNewSign_OldSign() const { return CurrentSign/OldSign; }

 void update_Sign();

 private:

 int CurrentSign, OldSign;

};


// little technical routine for move and measure construction
template <typename T1, typename T2>
inline std::string to_string (const T1 & x1, const T2 & x2) { 
 std::stringstream out;
 out<<x1<<"_"<<x2;
 return out.str();
}


#endif

