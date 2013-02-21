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
#include <triqs/utility/exceptions.hpp>
#include <triqs/arrays.hpp>

namespace triqs { namespace utility {

 namespace h5 = arrays::h5;

 template <typename T> struct type_to_number;
 template <> struct type_to_number<int> { static constexpr unsigned int value=1;};
 template <> struct type_to_number<long> { static constexpr unsigned int value=2;};
 template <> struct type_to_number<double> { static constexpr unsigned int value=3;};

 // ------------------------------------------------------

 struct _object {
  _object(){}
  _object(_object const &) = default;

  template<typename T> _object( T && x, std::string doc_ ): 
   p(std::make_shared<T>(std::forward<T>(x))),
   type_num(type_to_number<T>::value), 
   doc(doc_),
   h5_w ( [&x](h5::group_or_file const &F, std::string const &Name)->void { h5_write(F,Name, x);}) 
  {}
  
  
  std::shared_ptr<void> p;
  unsigned int type_num;
  std::string doc;
  std::function<void(h5::group_or_file const &, std::string const &)> h5_w;

  friend void h5_write ( h5::group_or_file F, std::string const & Name, _object const & obj){ obj.h5_w(F,Name); };

 };

 // ------------------------------------------------------
 // extraction 
 template<typename T> T extract(_object const &obj) {
  if (obj.type_num != type_to_number<T>::value) TRIQS_RUNTIME_ERROR<<"idiot !";return * static_cast<T*>(obj.p.get());
 }

 template<>
  double extract(_object const & obj) { 
   if (obj.type_num == type_to_number<double>::value) return * static_cast<double*>(obj.p.get());
   if (obj.type_num == type_to_number<int>::value) return extract<int>(obj);
   if (obj.type_num == type_to_number<long>::value) return extract<long>(obj);
   TRIQS_RUNTIME_ERROR<<"idiot !";
  }

 // ------------------------------------------------------
 class parameters {
  public : 

   parameters() = default;

   // put the default parameters ....
   //p(" Name",  DefaultValue, "Documentation")
   //
  //private: 
   typedef std::map<std::string, _object> map_t;
   map_t object_map, object_default_map;

   // template on a boolean ?? enable the = operator or delete ...
   template<bool IsConst>
    class ret_val {
     typedef typename std::conditional<IsConst, const parameters, parameters>::type param_t;

     param_t * const p;
     std::string const & key;
     template<typename T> T cast_impl() const { 
      auto it = p->object_map.find(key); 
      if (it!=p->object_map.end()) return extract<T>(it->second);
      auto it2 = p->object_default_map.find(key);
      if (it2==p->object_map.end()) TRIQS_RUNTIME_ERROR<<"Parameters : the key "<< key << " is not defined and not in the default";
      return extract<T>(it2->second);
     }
     public:
     ret_val (param_t * const p_, std::string const & key_ ): p(p_), key(key_) { }

#define CAST_OPERATOR(T) operator T () const{ return cast_impl<T>();}
     CAST_OPERATOR(long)
      CAST_OPERATOR(double)
      //CAST_OPERATOR(bool)
      //CAST_OPERATOR(std::string)
#undef CAST_OPERATOR

      template <typename RHS, bool C = IsConst> 
      typename std::enable_if<!C,  ret_val &>::type
      operator=(RHS && rhs) { 
       p->object_map[key] = _object(std::forward<RHS>(rhs),"");
       return *this; 
      }

     template <typename RHS, bool C = IsConst> 
      typename std::enable_if<C, ret_val &>::type
      operator=(RHS && rhs)  = delete;
    };

   friend class ret_val<true>; friend class ret_val<false>;

  public:

   ret_val<true> operator[](std::string const & k) const { return ret_val<true> (this,k); }
   ret_val<false> operator[](std::string const & k) { return ret_val<false> (this,k); }


 };

} } 
#endif

