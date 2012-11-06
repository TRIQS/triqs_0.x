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
#ifndef TRIQS_PYTHON_C_CONVERTERS_UNODERORED_MAP_H
#define TRIQS_PYTHON_C_CONVERTERS_UNODERORED_MAP_H
#include "../converters.hpp"

namespace triqs { namespace python_tools { 

 // boost::unordered_map<Key,Mapped> <----> python dict of conversion of Key, Mapped 
 template<typename Key, typename Mapped> 
  struct converter< boost::unordered_map<Key,Mapped> >  {

   typedef boost::unordered_map<Key,Mapped> C_type;
   typedef bpy::object P_type;

   static bool Py2C_is_possible (bpy::object d) { 
    try { 
     bpy::dict dic(d); bpy::list keys = dic.keys(), vals = dic.values();
     const size_t N = bpy::len(dic);
     for (size_t i=0; i<N; ++i) if ( (!converter<Key>::Py2C_is_possible(keys[i])) || (!converter<Mapped>::Py2C_is_possible(vals[i]))) throw false;
     return true;
    }
    catch(...) {return false;}
   }

   static C_type Py2C (bpy::object d) { 
    boost::unordered_map<Key,Mapped> res;
    bpy::dict dic(d);
    bpy::list keys = dic.keys(), vals = dic.values();
    const size_t N = bpy::len(dic);
    for (size_t i=0; i<N; ++i) {

     bool ok =  res.insert(std::make_pair(converter<Key>::Py2C(keys[i]), converter<Mapped>::Py2C(vals[i]))).second;
     if (!ok) TRIQS_RUNTIME_ERROR<<"This error should never happen !";
    }
    return res;
   }

   static P_type C2Py ( boost::unordered_map<Key,Mapped> const & m) {
    bpy::dict Res;
    for (typename boost::unordered_map<Key,Mapped>::const_iterator it = m.begin(); it !=m.end(); ++it) 
     Res[converter<Key>::C2Py (it->first)] = converter<Mapped>::C2Py (it->second);
    return Res;
   }
  };
}}
#endif
