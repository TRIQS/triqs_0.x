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
#ifndef TRIQS_UTILITY_PARAMS_H
#define TRIQS_UTILITY_PARAMS_H
#include <triqs/utility/first_include.hpp>
#include <string>
#include <complex>
#include <memory>
#include <map>
#include <boost/lexical_cast.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/serialization/map.hpp>
#include <triqs/utility/exceptions.hpp>
#include <triqs/utility/serialization.hpp>
#include <triqs/utility/typeid_name.hpp>
#include <triqs/arrays.hpp>

namespace triqs { namespace utility {

 using arrays::h5::group_or_file;
 using arrays::h5::write_attribute;
 using arrays::h5::h5_type_from_C;

 template<typename T> size_t type_hash_code () { return typeid(T).hash_code();}
 template<typename T> std::string type_name () { return demangle(typeid(T).name());}

 //template<typename T> std::string get_h5_type_name() { return type_name<T>();}

 // All the predefined cast of _object
#define TRIQS_UTIL_PARAM_PREDEFINED_CAST (int)(long)(double)(bool)(std::string)

 // 3 pts for hdf5 : read/write an attribute, get the type of data
 // Need for python : 
 // the constructor from PyObject * as in array
 // and to_ptyhon() methods, as in array ...

 // --------------- the opaque object ---------------------------------------
 // _object is a value : it makes deep copies in particular ...
 struct _object {

  std::shared_ptr<void> p;          // the object handled by the class
  size_t type_num;                  // id of the type, implementation dependent...
  std::string type_name_;           // for clear error message
  bool is_default_;
  std::string doc;                  // a documentation string, only used for default of parameters...
  std::function<_object()> clone_;  // clones the object : will be used to make copy of parameters !
  std::function<void(group_or_file const &, std::string const &)> h5_w; // the function to write in h5 !
  std::function<std::string()> serialize_; // for boost serialization ...
  std::function<void(std::ostream&)> print_;

  _object() {_init(); };

  template<typename T> _object( T * obj, std::string const & doc_ , bool is_default = false): 
   p(obj),
   type_num(type_hash_code<T>()), type_name_(type_name<T>()), 
   is_default_(is_default),
   doc(doc_),
   clone_( [obj,doc_]() { return _object( new T(*obj), doc_);}),
   h5_w ( [obj](group_or_file const &F, std::string const &Name)->void { h5_write(F,Name, *obj);}),
   //h5_w ( [obj](group_or_file const &F, std::string const &Name)->void { h5_write(F,Name, *obj); write_attribute(F,Name,get_h5_type_name<T>());}),
   serialize_([obj](){ return triqs::serialize(*obj);}) ,
   print_([obj](std::ostream & out ){out << *obj;})
   {_init();}

  _object(_object const & x ) { *this = x; _init(); } 
  _object(_object && c) { *this= std::move(c); _init();} 

  _object & operator = (_object && c) { swap(*this, c); return *this; }
  _object & operator = (_object const & x ) {*this = x.clone_(); return *this; }

  template <typename RHS> _object & operator=(RHS const & rhs) { *this =  _object (new RHS (rhs),""); return *this; }

  // special treatment for the const *, fall back to string  
  _object &  operator=(const char * rhs) { *this = std::string(rhs); return *this;}

  friend void swap (_object & a, _object & b) {  
#define SWAP(A)   swap(a.A,b.A)
   SWAP(p); std::SWAP(type_num); SWAP(type_name_); std::SWAP(is_default_);SWAP(doc); SWAP(clone_); SWAP(h5_w); SWAP(serialize_); SWAP(print_);
#undef SWAP
  }

  // implemented later, since need the extract function ...
#define CAST_OPERATOR(r, data, T) operator T () const;
  BOOST_PP_SEQ_FOR_EACH(CAST_OPERATOR, nil , TRIQS_UTIL_PARAM_PREDEFINED_CAST); 
#undef CAST_OPERATOR

  friend std::ostream & operator << (std::ostream & out, _object const & p) { p.print_(out); return out;}
  friend void h5_write ( group_or_file F, std::string const & Name, _object const & obj){ obj.h5_w(F,Name); };
  friend void h5_read ( group_or_file F, std::string const & Name, _object & obj);
  
  // ----- Boost serialisation 
  template<class Archive>
   void save(Archive & ar, const unsigned int version) const { 
    std::string s = serialize_();
    ar << boost::serialization::make_nvp("type_num", type_num);
    ar << boost::serialization::make_nvp("type_name", type_name_);
    ar << boost::serialization::make_nvp("seria_str", s);
   }
  template<class Archive>
   void load(Archive & ar, const unsigned int version) {
    std::string s, tnn; 
    size_t tn;
    ar >> boost::serialization::make_nvp("type_num", tn);
    ar >> boost::serialization::make_nvp("type_name", tnn);
    ar >> boost::serialization::make_nvp("seria_str", s);
    if (deserialization_fnts.find(tn)== deserialization_fnts.end())
     TRIQS_RUNTIME_ERROR << " Can not deserialize the type "<< tnn<< "\n Did you forget to register it ?";
    *this = deserialization_fnts[tn](s);
   }
  BOOST_SERIALIZATION_SPLIT_MEMBER();

  // ----- Unserialization (boost, h5) table : type_code -> deserialization function
  // init will register the most common types.
  static bool initialized;
  static void _init();
  static std::map<size_t, std::function<_object(std::string const &)>> deserialization_fnts;
  static std::map<size_t, std::function<_object(group_or_file const &, std::string const &)>> h5_read_fnts;
  static std::map<hid_t, std::pair<size_t,std::string> > h5_type_to_c_equivalence;

  template <typename T>
   static void register_type() {
    deserialization_fnts[type_hash_code<T>()] = [](std::string const &s) { return _object( new T( triqs::deserialize<T>(s)),"");};
    h5_read_fnts[type_hash_code<T>()] = [](group_or_file const &f,std::string const &s) ->_object { auto n = new T(); h5_read(f,s,*n); return _object(n,"");};
   }

  private: // native type
  template <typename T>
   static void register_native_type() {
    h5_type_to_c_equivalence[h5_type_from_C(T()).getId()] = std::make_pair(type_hash_code<T>(), type_name<T>());
    register_type<T>();
   }

 };// object class

 // --------------------- extraction  ---------------------------------
 template<typename T> T extract(_object const &obj);

 template<typename T> T lex_cast_from_string(_object const &obj) {
  std::string s = extract<std::string>(obj);
  try { return boost::lexical_cast<T>(s); } 
  catch(boost::bad_lexical_cast &) { TRIQS_RUNTIME_ERROR << " extraction : can not read the string "<<s <<" into a "<< type_name<T>(); }
 }

 template<typename T> T extract(_object const &obj) {
  // if T is not a string, and obj is string, attempt lexical_cast
  if ( (!std::is_same<T,std::string>::value) && (obj.type_num == type_hash_code<std::string>())) { return lex_cast_from_string<T>(obj); }
  if (obj.type_num != type_hash_code<T>()) 
   TRIQS_RUNTIME_ERROR<<"extraction : impossible : type mismatch. Got a "<<obj.type_name_<< " while I am supposed to extract a "<<type_name<T>();
  return * static_cast<T*>(obj.p.get());
 }

 template<> // specialize for double since we can make int -> conversion...
  double extract(_object const & obj) { 
   if (obj.type_num == type_hash_code<double>()) return * static_cast<double*>(obj.p.get());
   if (obj.type_num == type_hash_code<std::string>()) {return lex_cast_from_string<double>(obj); }
#define DECAYING_TYPE(T) if (obj.type_num == type_hash_code<T>()) return extract<T>(obj)
   DECAYING_TYPE(int);
   DECAYING_TYPE(unsigned int);
   DECAYING_TYPE(long);
   DECAYING_TYPE(unsigned long);
   DECAYING_TYPE(short);
   DECAYING_TYPE(unsigned short);
   DECAYING_TYPE(long long);
   DECAYING_TYPE(unsigned long long);
   DECAYING_TYPE(float);
#undef DECAYING_TYPE
   TRIQS_RUNTIME_ERROR<<"extraction : impossible : type mismatch. Got a "<<obj.type_name_<< " while I am supposed to extract a double";
  }

 // --------------- _object cast op implementation ---------------------------------------

#define CAST_OPERATOR(r, data, T) _object::operator T () const{ return extract<T>(*this);}
 BOOST_PP_SEQ_FOR_EACH(CAST_OPERATOR, nil , TRIQS_UTIL_PARAM_PREDEFINED_CAST); 
#undef CAST_OPERATOR

 // ------------------------------------------------------
 class parameters {
  public : 
   parameters() {};
   parameters (parameters const & other) = default;
   parameters (parameters && other) { swap(*this,other);}
   parameters & operator =  (parameters const & other)  = default;
   parameters & operator =  (parameters && other) { swap(*this,other); return *this;}
   friend void swap(parameters & a, parameters &b) { swap(a.object_map,b.object_map);} //swap(a.object_default_map,b.object_default_map);}

   // put the default parameters ....
   //p(" Name",  DefaultValue, "Documentation")
   //
 private: 
   typedef std::map<std::string, _object> map_t;
   map_t object_map; 

   friend class boost::serialization::access;
   template<class Archive>
    void serialize(Archive & ar, const unsigned int version) { ar & boost::serialization::make_nvp("object_map",object_map); }

 public:

   _object & operator[](std::string const & key) { return object_map[key];}

   _object const & operator[](std::string const & key) const { 
    auto it = object_map.find(key); 
    if ( it== object_map.end()) TRIQS_RUNTIME_ERROR<<"Key : "<< key<< " not found"; 
    return it->second;
   }

   friend void h5_write ( group_or_file F, std::string const & subgroup_name, parameters const & p){ 
    auto gr = F.create_group(subgroup_name);
    for (auto & pvp : p.object_map) h5_write(gr, pvp.first, pvp.second);
   }

   friend void h5_read ( group_or_file F, std::string const & subgroup_name, parameters & p);

   friend std::ostream & operator << (std::ostream & out, parameters const & p) { 
    out<< "{";
    for (auto & pvp : p.object_map) out<< pvp.first << " : " << pvp.second<< ", "; 
    return out<<"}";
   }
};

}} 
#endif
