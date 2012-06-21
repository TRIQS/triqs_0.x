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
#ifndef TRIQS_PYTHON_C_CONVERTERS_PAIR_H
#define TRIQS_PYTHON_C_CONVERTERS_PAIR_H
#include "../converters.hpp"

namespace triqs { namespace python_tools { 

 // std::pair<T1,T2> <----> python list or tuple of the conversion of T1,T2
 template<typename T1, typename T2> 
  struct converter<std::pair<T1,T2>, bpy::tuple> { 

   typedef std::pair<T1,T2> C_type;
   typedef bpy::tuple Py_type;

   static bool Py2C_is_possible (bpy::object l) { 
    try { 
     if (bpy::len(l)!=2) throw false;
     if (! (converter<T1>::Py2C_is_possible(l[0])) && (converter<T2>::Py2C_is_possible(l[1])) ) throw false;
     return true;
    }
    catch(...) {return false;}
   }

   static C_type Py2C (bpy::object l) { 
    bpy::list L(l);
    const size_t N = bpy::len(L);
    if (N!=2) TRIQS_RUNTIME_ERROR<<"Object"<<object_to_string(l)<<" can not be converted to a std::pair : lengths mismatch "; 
    return std::make_pair( converter<T1>::invoke(l[0]), converter<T2>::invoke(l[1]) );
   }

   static bpy::object C2Py (C_type const & p) { return bpy::make_tuple( converter<T1>::invoke (p.first), converter<T2>::invoke (p.second));}
  };
}}
#endif
