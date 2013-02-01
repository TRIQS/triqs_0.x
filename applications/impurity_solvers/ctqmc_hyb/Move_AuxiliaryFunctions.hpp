
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

#ifndef MOVE_AUXILIARY_H_ehfkh
#define MOVE_AUXILIARY_H_ehfkh

#include "Configuration.hpp"

const double EPSILON = 1.e-13;

#ifdef DEBUG

inline void print_det(Configuration::DET_TYPE * det){
   cout<<"*****************"<<endl;
   cout<<"Det: NumC = "<<det->size()<<endl;
   if (det->size()==0) {cout<<"VIDE\n*****************"<<endl; return;}
   for (Configuration::DET_TYPE::Cdagger_iterator p= det->Cdagger_begin(); (p != det->Cdagger_end()) ; ++p) 
     cout<<" C_dag "<< p->tau<<"  "<<p->Op->name<< endl;
   for (Configuration::DET_TYPE::C_iterator p= det->C_begin(); (p != det->C_end())  ; ++p) 
     cout<<" C "<< p->tau<<"  "<<p->Op->name<< endl;
   cout<<"*****************"<<endl;
} 

inline void print_det_reverse(Configuration::DET_TYPE * det){
   cout<<"*****************"<<endl;
   cout<<"Det: NumC = "<<det->size()<<endl;
   if (det->size()==0) {cout<<"VIDE\n*****************"<<endl; return;}
   for (Configuration::DET_TYPE::Cdagger_reverse_iterator p= det->Cdagger_rbegin(); (p != det->Cdagger_rend()) ; ++p) 
     cout<<" C_dag "<< p->tau<<"  "<<p->Op->name<< endl;
   for (Configuration::DET_TYPE::C_reverse_iterator p= det->C_rbegin(); (p != det->C_rend())  ; ++p) 
     cout<<" C "<< p->tau<<"  "<<p->Op->name<< endl;
   cout<<"*****************"<<endl;
} 
#endif
 

#endif
