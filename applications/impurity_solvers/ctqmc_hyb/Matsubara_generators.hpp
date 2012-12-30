
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

#ifndef MATSU_GEN_H
#define MATSU_GEN_H

#include "OP_Tools/all.hpp"

/** 
    Generates the Exp(i(2n+1) Pi/Beta)
*/
class Matsubara_Expo_Generator
{
public:
  /**
   */
  Matsubara_Expo_Generator(double Beta, double t,COMPLEX coef=1) : 
    Start(exp( I * Pi/Beta * t)),Step(Start*Start),Coef(coef){reset();}
  
  /**
   */
  inline void reset(){ om = Coef*Start;} 
  
  /**
   */
  inline COMPLEX operator()(){return om;}
  
  /**
     Prefix notation
   */
  inline void operator++() { om *= Step;}
    
protected : 
  const COMPLEX Start,Step,Coef;
  COMPLEX om;
};


/** 
    Generates the i(2n+1) Pi/Beta
*/
class Matsubara_Omega_Generator
{
public:
  /**
   */
  Matsubara_Omega_Generator(double Beta, COMPLEX shift=0) : 
    Start( I * Pi/Beta + shift),Step(2*(Start-shift)){reset();}
  
  /**
   */
  inline void reset(){ om = Start;} 
  
  /**
   */
  inline COMPLEX operator()(){return om;}
  
  /**
   */
  inline void operator++() { om += Step;}
  
protected : 
  COMPLEX Start;
  COMPLEX Step,om;
};

#endif

