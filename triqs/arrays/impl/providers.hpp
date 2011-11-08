
/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011 by O. Parcollet
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

#ifndef TRIQS_ARRAYS_PROVIDERS_H
#define TRIQS_ARRAYS_PROVIDERS_H
#include "./common.hpp"
#include "./compound_assign.hpp"

namespace triqs { namespace arrays {  namespace providers { 

 // add infix operation with CRTP
 template<typename Derived> struct compound_assign_ops { 

  template<typename RHS> 
   Derived & operator +=(RHS const & rhs) { 
    Derived & self = static_cast<Derived & >(*this); 
    compound_assignment<Derived,RHS,'A'>(self,rhs);
    return self;  
   }

  template<typename RHS> 
   Derived & operator -=(RHS const & rhs) { 
    Derived & self = static_cast<Derived & >(*this); 
    compound_assignment<Derived,RHS,'S'>(self,rhs);
    return self;  
   }

  template<typename RHS> 
   Derived & operator *=(RHS const & rhs) { 
    Derived & self = static_cast<Derived & >(*this); 
    compound_assignment<Derived,RHS,'M'>(self,rhs);
    return self;  
   }

  template<typename RHS> 
   Derived & operator /=(RHS const & rhs) { 
    Derived & self = static_cast<Derived & >(*this); 
    compound_assignment<Derived,RHS,'D'>(self,rhs);
    return self;  
   }

 };
}}} //namespace triqs::arrays::providers
#endif

