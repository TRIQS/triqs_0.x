/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013 by H. Hafermann, O. Parcollet
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
#ifndef TRIQS_UTILITY_PARAMS_DEFAULT_H
#define TRIQS_UTILITY_PARAMS_DEFAULT_H
#include "./opaque_object_h5.hpp"
namespace triqs { namespace utility {
 /**
  * DOC TO BE WRITTEN
  */
 class parameters_defaults {
  public :
   parameters_defaults() {};
   parameters_defaults (parameters_defaults const & other) = default;
   parameters_defaults (parameters_defaults && other) { swap(*this,other);}
   parameters_defaults & operator =  (parameters_defaults const & other)  = default;
   parameters_defaults & operator =  (parameters_defaults && other) { swap(*this,other); return *this;}
   friend void swap(parameters_defaults & a, parameters_defaults &b) 
   { swap(a.object_map,b.object_map); swap(a.documentation, b.documentation);swap(a.is_optional, b.is_optional);  }

  private:
   typedef std::map<std::string, _object> map_t;
   map_t object_map;
   std::map<std::string, std::string> documentation;
   std::map<std::string, bool > is_optional;

   friend class boost::serialization::access;
   template<class Archive>
    void serialize(Archive & ar, const unsigned int version) { ar & boost::serialization::make_nvp("object_map",object_map); }

   struct _inserter {
    parameters_defaults * p; bool opt;
    _inserter(parameters_defaults *p_, bool opt_) : p(p_), opt(opt_) {}
    template<typename T> _inserter operator()(std::string const & key, T && def_val, std::string const & doc) {
     p->object_map[key] = std::forward<T>(def_val);
     p->documentation[key] = doc;
     p->is_optional[key] = opt;
     return *this;
    }

   };
   friend struct _inserter;

   template<typename T> const T getter(std::map<std::string,T> const & m, std::string const & key) const {
    auto it = m.find(key); assert(it !=m.end()); return it->second;
   }

  public:

   typedef map_t::const_iterator const_iterator;
   typedef map_t::iterator iterator;
   const_iterator begin() const { return object_map.begin();}
   const_iterator end()  const { return object_map.end();}
   iterator begin() { return object_map.begin();}
   iterator end() { return object_map.end();}

   bool has_key(std::string const & key) const { return object_map.find(key) != object_map.end();}
   bool is_required(std::string const & key) const { return (has_key(key) && (! getter(this->is_optional,key)));}
   std::string doc(std::string const & key) const { return (has_key(key) ? getter(this->documentation,key) : "");}

   template<typename T>
    _inserter optional (std::string const & key, T && def_val, std::string const & doc) {
     return _inserter(this, true)(key,std::forward<T>(def_val), doc);
    }

   template<typename T>
    _inserter required (std::string const & key, T && def_val, std::string const & doc) {
     return _inserter(this, false)(key,std::forward<T>(def_val), doc);
    }

   ///
   _object const & operator[](std::string const & key) const {
    auto it = object_map.find(key);
    if ( it== object_map.end()) TRIQS_RUNTIME_ERROR<<"Key : "<< key<< " not found";
    return it->second;
   }

   friend std::ostream & operator << (std::ostream & out, parameters_defaults const & p) {
    out<< "{";
    for (auto & pvp : p.object_map) out<< pvp.first << " : " << pvp.second<< ", ";
    return out<<"}";
   }

 };

}}
#endif
