
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

#ifndef TOOLS_H
#define TOOLS_H

#include <string>
#include <sstream>

// DEBUG TOOL ....

#ifdef __GNUC__
#include <cxxabi.h> // gcc only ?
#endif

namespace triqs { namespace arrays { namespace Tools { 

 template<typename T>
  std::string typeid_name(T const & A) { 
   std::stringstream fs;
#ifdef __GNUC__
   int status;
   fs<<abi::__cxa_demangle(typeid(A).name(), 0, 0, &status);
#else 
   fs<<typeid(A).name();
#endif
   return fs.str();
  }

};};
#endif
