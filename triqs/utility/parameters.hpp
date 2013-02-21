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

 template<typename T> size_t type_hash_code () { return typeid(T).hash_code();}
 template<typename T> std::string type_name () { return demangle(typeid(T).name());}

 template<typename T> std::string get_h5_type_name() { return type_name<T>();}

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
  std::string doc;                  // a documentation string, only used for default of parameters...
  std::function<_object()> clone_;  // clones the object : will be used to make copy of parameters !
  std::function<void(group_or_file const &, std::string const &)> h5_w; // the function to write in h5 !
  std::function<std::string()> serialize_; // for boost serialization ...
  std::function<void(std::ostream&)> print_;

  _object() {};

  template<typename T> _object( T * obj, std::string const & doc_ ): 
   p(obj),
   type_num(type_hash_code<T>()), type_name_(type_name<T>()), 
   doc(doc_),
   clone_( [obj,doc_]() { return _object( new T(*obj), doc_);}),
   h5_w ( [obj](group_or_file const &F, std::string const &Name)->void { h5_write(F,Name, *obj);}),
   //h5_w ( [obj](group_or_file const &F, std::string const &Name)->void { h5_write(F,Name, *obj); write_attribute(F,Name,get_h5_type_name<T>());}),
   serialize_([obj](){ return triqs::serialize(*obj);}) ,
   print_([obj](std::ostream & out ){out << *obj;})
   {}

  friend std::ostream & operator << (std::ostream & out, _object const & p) { p.print_(out); return out;}
    
  friend void h5_write ( group_or_file F, std::string const & Name, _object const & obj){ obj.h5_w(F,Name); };
 
  friend void h5_read ( group_or_file F, std::string const & Name, _object & obj){ 
   //to try only. Need to make h5_code function for all object....
   // write, read an attribute with the type ? 
   obj = h5_read_fnts[type_hash_code<double>()](F,Name);
  }

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

  // table :  code -> _object ( string) 
  static std::map<size_t, std::function<_object(std::string const &)>> deserialization_fnts;
  static std::map<size_t, std::function<_object(group_or_file const &, std::string const &)>> h5_read_fnts;

  template <typename T>
   static void register_type() {
    deserialization_fnts[type_hash_code<T>()] = [](std::string const &s) { return _object( new T( triqs::deserialize<T>(s)),"");};
     h5_read_fnts[type_hash_code<T>()] = [](group_or_file const &f,std::string const &s) ->_object { auto n = new T(); h5_read(f,s,*n); return _object(n,"");};
   }
 };


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

 // ------------------------------------------------------
 class parameters {
  public : 

   parameters() { _init(); };

   parameters (parameters const & other) {
    for (auto & pvp : other.object_map)         { object_map.insert( std::make_pair(pvp.first,pvp.second.clone_()));} 
    for (auto & pvp : other.object_default_map) { object_default_map.insert( std::make_pair(pvp.first,pvp.second.clone_()));} 
    _init();
   }

   parameters (parameters && other) { swap(*this,other);}

   friend void swap(parameters & a, parameters &b) { swap(a.object_map,b.object_map);swap(a.object_default_map,b.object_default_map);}

   parameters & operator =  (parameters const & other)  = default;
   parameters & operator =  (parameters && other) { swap(*this,other); return *this;}

   // put the default parameters ....
   //p(" Name",  DefaultValue, "Documentation")
   //
  private: 
   typedef std::map<std::string, _object> map_t;
   map_t object_map, object_default_map;

   // return class of [] operator
   template<bool IsConst> class ret_val {
    typedef typename std::conditional<IsConst, const parameters, parameters>::type param_t;
    param_t * const p;
    std::string const & key;
    _object get_object() const { 
     auto it = p->object_map.find(key); 
     if (it!=p->object_map.end()) return it->second;
     auto it2 = p->object_default_map.find(key);
     if (it2==p->object_map.end()) TRIQS_RUNTIME_ERROR<<"Parameters : the key "<< key << " is not defined and not in the default";
     return it2->second;
    }
    public:
    ret_val (param_t * const p_, std::string const & key_ ): p(p_), key(key_) { }

#define ALLOWED_CAST (int)(long)(double)(bool)(std::string)
#define CAST_OPERATOR(r, data, T) operator T () const{ return extract<T>(get_object());}
    BOOST_PP_SEQ_FOR_EACH(CAST_OPERATOR, nil , ALLOWED_CAST); 
#undef CAST_OPERATOR

    template <typename RHS, bool C = IsConst> 
     typename std::enable_if<!C,  ret_val &>::type
     operator=(RHS const & rhs) {
      // if (p->object_map.find(key) != p->object_map.end())// do I need this ? 
      p->object_map.erase(key);
      p->object_map.insert(std::make_pair(key, _object (new RHS (rhs),"")));
      return *this; 
     }

    template <bool C = IsConst> // special treatment for the const *  
     typename std::enable_if<!C,  ret_val &>::type
     operator=(const char * rhs) { (*this) = std::string(rhs); return *this;}

    friend std::ostream & operator << (std::ostream & out, ret_val const & x) { return out<< x.get_object();}
   };

   friend class ret_val<true>; friend class ret_val<false>;

   friend class boost::serialization::access;
   template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
     ar & boost::serialization::make_nvp("object_map",object_map);
     ar & boost::serialization::make_nvp("object_default_map",object_default_map);
    }

   static bool initialized;
   static void _init() { 
    if (initialized) return;
#define REGISTER_UNSERIALIZE_FUN(r, data, T) _object::register_type<T>(); 
    BOOST_PP_SEQ_FOR_EACH(REGISTER_UNSERIALIZE_FUN, nil , ALLOWED_CAST); 
#undef REGISTER_UNSERIALIZE_FUN
    _object::register_type<parameters>(); 
    initialized = true;
   }

  public:

   ret_val<true> operator[](std::string const & k) const { return ret_val<true> (this,k); }
   ret_val<false> operator[](std::string const & k) { return ret_val<false> (this,k); }

   friend void h5_write ( group_or_file F, std::string const & subgroup_name, parameters const & p){ 
    auto gr = F.create_group(subgroup_name);
    for (auto & pvp : p.object_map) h5_write(gr, pvp.first, pvp.second);
   }

   friend void h5_read ( group_or_file F, std::string const & subgroup_name, parameters & p){ 
    auto gr = F.open_group(subgroup_name);
    _object obj; 
    // need to loop over subgroup
    h5_read(gr,"a",obj); p["a"] = obj;
    //for (auto & pvp : p.object_map) h5_read(gr, pvp.first, pvp.second);
   }

   friend std::ostream & operator << (std::ostream & out, parameters const & p) { 
    out<< "{";
    for (auto & pvp : p.object_map) out<< pvp.first << " : " << pvp.second<< ", "; //<< std::endl ;
    return out<<"}";
   }
 };

}} 
#endif
