
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

#ifndef DETMANIP_TESTER_H
#define DETMANIP_TESTER_H

#include "./det_manip.hpp"
namespace OldVersion { 
#include "./detManip.old.hpp"
}

#define DET_MANIP_CHECKER1
#define DET_MANIP_CHECKER

template<class DELTATYPE,class VALTYPE,class TAUTYPE_PTR>
struct FunctionAdapt { 

 typedef TAUTYPE_PTR argument_type;
 typedef VALTYPE return_type;
 DELTATYPE delta;
 FunctionAdapt(DELTATYPE const & delta_): delta(delta_){}
 return_type operator()(argument_type const & x, argument_type const & y) { return delta(x,y);}
};

template<class DELTATYPE,class VALTYPE,class TAUTYPE_PTR>
class detManip : public OldVersion::detManip<DELTATYPE,VALTYPE,TAUTYPE_PTR> 
{
 typedef FunctionAdapt<DELTATYPE,VALTYPE,TAUTYPE_PTR> Function;
 det_manip<Function> d2;
 OldVersion::detManip<DELTATYPE,VALTYPE,TAUTYPE_PTR> d1;
 static const bool use_d1 = true;
 static const bool use_d2 = true;

 //double res,res2;
 typedef OldVersion::detManip<DELTATYPE,VALTYPE,TAUTYPE_PTR> BaseType;

 public:

 detManip(const DELTATYPE & Delta,int Nmax_,VALTYPE Eta_): 
  BaseType(Delta, Nmax_, Eta_),d1(Delta,Nmax_,Eta_), 
  d2(Function(Delta),Nmax_) {}

 VALTYPE try_insert(int i0, int j0, TAUTYPE_PTR tau, TAUTYPE_PTR taup) {
  double res =  BaseType::try_insert(i0,j0,tau,taup);
#ifdef DET_MANIP_CHECKER
  double res2 =  d2.try_insert(i0-1,j0-1,tau,taup);
  CHECK_OR_THROW((abs(res-res2)<1.e-12), "Error in insert"<<i0<<j0<<res<<" "<<res2); //<<tau<<taup);
#endif
#ifdef DET_MANIP_CHECKER1
  d1.try_insert(i0,j0,tau,taup);
#endif
  return res;
 }

 VALTYPE try_remove(int i0, int j0) {
 double res =  BaseType::try_remove(i0,j0);
#ifdef DET_MANIP_CHECKER
  double res2 =  d2.try_remove(i0-1,j0-1);
  CHECK_OR_THROW((abs(res-res2)<1.e-12), "Error in remove"<<i0<<j0);
#endif
#ifdef DET_MANIP_CHECKER1
  d1.try_remove(i0,j0);
#endif
  return res;
 }

 VALTYPE try_change_col(int i0, TAUTYPE_PTR Taup) {
  double res =  BaseType::try_change_col(i0,Taup);
#ifdef DET_MANIP_CHECKER
  double res2 =  d2.try_change_col(i0-1,Taup);
    CHECK_OR_THROW((abs(res-res2)<1.e-12), "Error in change_col"<<j0);//<<Taup);
#endif
#ifdef DET_MANIP_CHECKER1
  d1.try_change_col(i0,Taup);
#endif
    return res;
 }

 VALTYPE try_change_row(int i0, TAUTYPE_PTR Tau) {
  double res =  BaseType::try_change_row(i0,Tau);
#ifdef DET_MANIP_CHECKER
  double res2 =  d2.try_change_row(i0-1,Tau);
  CHECK_OR_THROW((abs(res-res2)<1.e-12), "Error in change_row"<<i0); // <<Tau);
#endif
#ifdef DET_MANIP_CHECKER1
  d1.try_change_row(i0,Tau);
#endif
  return res;
 }

 void accept_move(){
  BaseType::accept_move();
#ifdef DET_MANIP_CHECKER
  d2.complete_operation();
#endif
#ifdef DET_MANIP_CHECKER1
  d1.accept_move();
#endif
 }

 inline int roll_matrix(typename BaseType::RollDirection roll) { 
 int res =  BaseType::roll_matrix(roll);
#ifdef DET_MANIP_CHECKER
 int res2 =  d2.roll_matrix(typename det_manip<Function>::RollDirection(roll));
  CHECK_OR_THROW((abs(res-res2)<1.e-12), "Error in roll");
#endif
#ifdef DET_MANIP_CHECKER1
 d1.roll_matrix(roll);
#endif
  return res;
 }

 inline VALTYPE determinant() const {
  double res = BaseType::determinant();
#ifdef DET_MANIP_CHECKER
  double res2 = d2.determinant();
  CHECK_OR_THROW((abs(res-res2)<1.e-12), "Error in det");
#endif
  return res;
 }
};

#endif


