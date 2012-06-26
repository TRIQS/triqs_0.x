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
#ifndef TRIQS_PYTHON_C_CONVERTERS_VECTOR_H 
#define TRIQS_PYTHON_C_CONVERTERS_VECTOR_H
#include "../converters.hpp"

namespace triqs { namespace python_tools { 

 // std::vector<T> <----> python list or tuple of the conversion of T
 template<typename T> struct converter<std::vector<T> > { 

  typedef std::vector<T> C_type;
  typedef bpy::list P_type;

  static bool Py2C_is_possible (bpy::object l) { 
   try { 
    P_type L(l); // check it is a list...
    const size_t N = bpy::len(L);
    for (size_t i=0; i<N; ++i) if (!converter<T>::Py2C_is_possible(l[i])) throw false;
    return true;
   }
   catch(...) {return false;}
  }

  static C_type Py2C (bpy::object l) { 
   P_type L(l);
   const size_t N = bpy::len(L);
   std::vector<T> res; res.reserve(N);
   for (size_t i=0; i<N; ++i) res.push_back( converter<T>::Py2C(l[i]));
   return res;
  }

  static P_type C2Py (C_type const &v) {
   P_type L;
   for (size_t i=0; i<v.size(); ++i) L.append( converter<T>::C2Py(v[i])); 
   return L;
  }
 };

}}
#endif
