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

 template<typename T> constexpr size_t type_hash_code () { return typeid(T).hash_code();}
 template<typename T> std::string type_name () { return demangle(typeid(T).name());}

 // All the predefined cast of _object
#define TRIQS_UTIL_PARAM_PREDEFINED_CAST (int)(long)(double)(bool)(std::string)

 // --------------- the opaque object ---------------------------------------
 // _object is a value : it makes deep copies in particular ...
 struct _object {

  std::shared_ptr<void> p;          // the object handled by the class
  size_t type_num;                  // id of the type, implementation dependent...
  std::string type_name_;           // for clear error message
  std::function<_object()> clone_;  // clones the object : will be used to make copy of parameters !
  std::function<void(h5::group const &, std::string const &)> h5_w; // the function to write in h5 !
  std::function<std::string()> serialize_; // for boost serialization ...
  std::function<void(std::ostream&)> print_;

  _object() {_init(); };

  // pass a const & or a && (in which case a move will be used, provided that T has a move constructor
  template<typename T> _object( T && obj, const char *){ // const char *  is there to disambigue with copy construction
   delegate_constructor( new typename std::remove_reference<T>::type(std::forward<T>(obj)));
  }

  private:
  template<typename T> void delegate_constructor( T * obj) {
   p = std::shared_ptr<void> (obj);
   type_num = type_hash_code<T>();
   type_name_ = type_name<T>();
   clone_ =  [obj]() { return _object( *obj, "");} ;
   h5_w =  [obj](h5::group const &F, std::string const &Name)->void { h5_write(F,Name, *obj);};
   serialize_ = [obj](){ return triqs::serialize(*obj);};
   print_ = [obj](std::ostream & out ){out << *obj;};
   _init();
  }

  public:

  _object(_object const & x ) { *this = x; _init(); }
  _object(_object && c) { *this= std::move(c); _init();}

  friend void swap (_object & a, _object & b) {
#define SWAP(A)   swap(a.A,b.A)
   SWAP(p); std::SWAP(type_num); SWAP(type_name_); SWAP(clone_); SWAP(h5_w); SWAP(serialize_); SWAP(print_);
#undef SWAP
  }

  _object & operator = (_object && c) { swap(*this, c); return *this; }
  _object & operator = (_object const & x ) {*this = x.clone_(); return *this; }

  template <typename RHS> _object & operator=(RHS && rhs) {
   *this =  _object (std::forward<RHS>(rhs),""); 
   register_type<typename std::remove_reference<RHS>::type>::invoke(); 
   return *this; 
  }

  // special treatment for the const *, fall back to string
  _object &  operator=(const char * rhs) { *this = std::string(rhs); return *this;}

  // implemented later, since need the extract function ...
#define CAST_OPERATOR(r, data, T) operator T () const;
  BOOST_PP_SEQ_FOR_EACH(CAST_OPERATOR, nil , TRIQS_UTIL_PARAM_PREDEFINED_CAST);
#undef CAST_OPERATOR

  friend std::ostream & operator << (std::ostream & out, _object const & p) { p.print_(out); return out;}
  friend void h5_write ( h5::group F, std::string const & Name, _object const & obj){ obj.h5_w(F,Name); };
  friend void h5_read ( h5::group F, std::string const & Name, _object & obj);

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
  static std::map<size_t, std::function<_object(h5::group const &, std::string const &)>> h5_read_fnts;
  static std::map<hid_t, std::pair<size_t,std::string> > h5_type_to_c_equivalence;
  static std::map<std::pair<size_t,int>, size_t > hash_code_element_rank_2_hash_code_array;
  static std::map<std::string, size_t > h5_scheme_to_code;

  template <typename T> struct register_type;

  private: // native type
  template <typename T>
   static void register_native_type() {
    if (register_type<T>::invoke()) return;
    h5_type_to_c_equivalence[h5::h5_type_from_C(T()).getId()] = std::make_pair(type_hash_code<T>(), type_name<T>());
   }

 };// object class

 // --------------------- type registering  ---------------------------------

 template <typename T> struct _object::register_type {
  static bool invoke() { // returns true if already registered
   size_t code = type_hash_code<T>();
   if (deserialization_fnts.find(code)!= deserialization_fnts.end()) return true;
   deserialization_fnts[code] = [](std::string const &s) { return _object( triqs::deserialize<T>(s),"");};
   h5_read_fnts[code] = [](h5::group const &f,std::string const &s) ->_object { T n; h5_read(f,s,n); return _object(std::move(n),"");};
   auto h5_scheme = get_triqs_hdf5_data_scheme(T());
   if (h5_scheme != "") h5_scheme_to_code[h5_scheme] = code;
   std::cerr  << " registering " << type_name<T>() << "h5 scehme "<< h5_scheme<< std::endl ;
   return false;
  }
 };

 // special case for array, we need to fill one more table
 template<typename T, int R> struct _object::register_type<arrays::array<T,R>>{ 
  static bool invoke() {
   typedef arrays::array<T,R> A;
   if (deserialization_fnts.find(type_hash_code<A>())!= deserialization_fnts.end()) return true;
   //std::cerr  << " registering " << type_name<A>() << std::endl ;
   deserialization_fnts[type_hash_code<A>()] = [](std::string const &s) { return _object( triqs::deserialize<A>(s),"");};
   h5_read_fnts[type_hash_code<A>()] = [](h5::group const &f,std::string const &s) ->_object { A n; h5_read(f,s,n); return _object(std::move(n),"");};
   hash_code_element_rank_2_hash_code_array[std::make_pair(type_hash_code<T>(), R)] = type_hash_code<A>();
   return false;
  }
 };

 // --------------------- arithmetic operations are deleted for _object  ---------------------------------

#define DELETE_OP(op)\
 template <typename LHS, typename RHS> \
 typename std::enable_if< std::is_same<LHS,_object>::value || std::is_same<RHS,_object>::value>::type\
 operator op( LHS const &, RHS const &) = delete;
 DELETE_OP(+); DELETE_OP(-); DELETE_OP(*); DELETE_OP(/);
#undef DELETE_OP

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
   //DECAYING_TYPE(unsigned int);
   DECAYING_TYPE(long);
   //DECAYING_TYPE(unsigned long);
   DECAYING_TYPE(short);
   //DECAYING_TYPE(unsigned short);
   DECAYING_TYPE(long long);
   //DECAYING_TYPE(unsigned long long);
   //DECAYING_TYPE(float);
#undef DECAYING_TYPE
   TRIQS_RUNTIME_ERROR<<"extraction : impossible : type mismatch. Got a "<<obj.type_name_<< " while I am supposed to extract a double";
  }

 // --------------- _object cast op implementation ---------------------------------------

#define CAST_OPERATOR(r, data, T) _object::operator T () const{ return extract<T>(*this);}
 BOOST_PP_SEQ_FOR_EACH(CAST_OPERATOR, nil , TRIQS_UTIL_PARAM_PREDEFINED_CAST);
#undef CAST_OPERATOR

 // ------------------------------------------------------

 /**
  * DOC TO BE WRITTEN
  */
 class parameters {
  public :
   parameters() {};
   parameters (parameters const & other) = default;
   parameters (parameters && other) { swap(*this,other);}
   parameters & operator =  (parameters const & other)  = default;
   parameters & operator =  (parameters && other) { swap(*this,other); return *this;}
   friend void swap(parameters & a, parameters &b) { swap(a.object_map,b.object_map);} 

  private:
   typedef std::map<std::string, _object> map_t;
   map_t object_map;
   std::map<std::string, std::string> documentation;

   friend class boost::serialization::access;
   template<class Archive>
    void serialize(Archive & ar, const unsigned int version) { ar & boost::serialization::make_nvp("object_map",object_map); }

   struct _inserter { 
    parameters * p; 
    _inserter(parameters *p_) : p(p_) {}
    template<typename T> _inserter operator()(std::string const & key, T && def_val, std::string const & doc) {
      p->object_map[key] = std::forward<T>(def_val); 
      p->documentation[key] = doc;
      return *this;
     }
   };
   friend struct _inserter;

  public:

   typedef map_t::const_iterator const_iterator;
   typedef map_t::iterator iterator;
   const_iterator begin() const { return object_map.begin();}
   const_iterator end()  const { return object_map.end();}
   iterator begin() { return object_map.begin();}
   iterator end() { return object_map.end();}

   bool has_key(std::string const & k) const { return object_map.find(k) != object_map.end();}

   /// Insert p[key] = def_val and keep the corresponding documentation ...
   template<typename T> 
    _inserter insert(std::string const & key, T && def_val, std::string const & doc) { 
     return _inserter(this)(key,std::forward<T>(def_val), doc);
    }

   ///
   _object & operator[](std::string const & key) { return object_map[key];}

   ///
   _object const & operator[](std::string const & key) const {
    auto it = object_map.find(key);
    if ( it== object_map.end()) TRIQS_RUNTIME_ERROR<<"Key : "<< key<< " not found";
    return it->second;
   }

   friend void h5_write ( h5::group F, std::string const & subgroup_name, parameters const & p){
    auto gr = F.create_group(subgroup_name);
    for (auto & pvp : p.object_map) h5_write(gr, pvp.first, pvp.second);
   }

   friend void h5_read ( h5::group F, std::string const & subgroup_name, parameters & p);

   friend std::ostream & operator << (std::ostream & out, parameters const & p) {
    out<< "{";
    for (auto & pvp : p.object_map) out<< pvp.first << " : " << pvp.second<< ", ";
    return out<<"}";
   }

   /**
    * Register a type for conversion, serialization and h5 read/write.
    * Note : can be called multiple times (no effect for second and later call).
    * Note : this is automatically called when putting an object in parameters
    */
   template<typename T> static void register_type() { _object::register_type<T>::invoke();}

   /** 
    * Update all keys with pdef, possibly overwriting without any check (python like behaviour).
    */
   void update(parameters const & pdef);

   /// Flags controlling the update_with_default function
   //static constexpr ull_t strict_type_check = 1ull;              // Type check is strict. alsway true now
   static constexpr ull_t reject_key_without_default = 1ull<<2;

   /** 
    * Update with a default parameters set.
    * If a key is present in pdef and not in this, add the pdef value to this
    * If a key is present in both, do no change it in this, but check that type are the same (to be smoothed)
    * If a key is present in this, and not in pdef, and reject_key_without_default is passed, then raise triqs::runtime_error exception
    */
   void update_with_defaults(parameters const & pdef, ull_t flag =0);

 };

}}
#endif
