
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

#include "./stack_trace.hpp"
#include <sstream>

#ifndef TRIQS_TRACE_MAX_FRAMES
#define TRIQS_TRACE_MAX_FRAMES 50
#endif

#ifdef __GNUC__

#include <execinfo.h>
#include <cxxabi.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include "./typeid_name.hpp"

namespace triqs { namespace utility {  

 std::string stack_trace() {
  std::ostringstream buffer;

  void * stack[TRIQS_TRACE_MAX_FRAMES + 1];
  std::size_t depth = backtrace(stack, TRIQS_TRACE_MAX_FRAMES + 1);
  if (!depth)
   buffer << "  empty  " << std::endl;
  else {
   char * * symbols = backtrace_symbols(stack, depth);
   for (std::size_t i = 1; i < depth; ++i) {
    std::string symbol = symbols[i];
    std::vector<std::string> strs;
    boost::split(strs, symbol, boost::is_any_of("\t ()+"), boost::token_compress_on);
    for (std::vector<std::string>::const_iterator it = strs.begin(); it !=strs.end(); ++it)
     buffer << " "<<triqs::utility::demangle(*it);
    buffer << std::endl;
    //buffer << ": " << symbol << std::endl;
   }
   free(symbols);
  }
  return buffer.str();
 }
}}
#else

namespace triqs { namespace utility {  

 std::string stack_trace() { return std::string("stacktrace only available in gcc");}

}}

#endif


