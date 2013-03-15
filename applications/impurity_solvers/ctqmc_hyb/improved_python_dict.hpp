
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

#ifndef IMPROVED_PYTHON_DICT_H
#define IMPROVED_PYTHON_DICT_H
#include <boost/python.hpp>

namespace triqs { namespace python_tools {

 /**
  * A dressed version of the python::dict, which removes the need of extract for simple types.
  * To avoid ambiguous situations the extraction to dicts or lists must still be done explicitly
  * using the direct access to the underlying dict with dict()
  *
  * Examples:
  * 
  * int i = D['a'];                          // same as int i = extract<int>(D.dict()['a'])
  * dict d = extract<dict>(D.dict()['b']);   // d is a view of the dict D['b']
  * list l = extract<list>(D.dict()['c']);   // l is a view of the list D['c']
  * 
  * Details:
  *
  * dict() : returns the underlying python::dict
  * D['aa'] is a proxy that can be casted to anything (it applies extract).
  * D['aa'] = x modifies the dictionnary. 
  * If 'aa' is not in the dictionary, an exception is thrown.
  *
  * D.value_or_default('aaa',default) : returns the value or the default value if the value is not in the dict
  *
  * This class has the same interface as alps::mcparams for alps compatibility. 
  */
 class improved_python_dict {  
  boost::python::dict _d;
  class val; class cval; // forward, to keep the public part first... just convenience

  public : 
  
  /// d must be a dict. No copy is made of d, the class acts on the dictionnary
  improved_python_dict (boost::python::object d) {
    boost::python::extract<boost::python::dict> ex(d);
    if (!ex.check()) throw "improved_python_dict : construction : the object is not a dict !";
    _d = ex();
  };

  val operator[](std::string const & k) { defined_or_throw(k); return val (_d,k); }
  cval operator[](std::string const & k) const { defined_or_throw(k); return cval (_d,k); }

  boost::python::object dict() const { return _d;}

  template<class T>
  T value_or_default(std::string const & k, T const & x) const { return (defined(k) ? boost::python::extract<T>(_d[k]) : x); }

  /// const char * is transformed into a string
  std::string value_or_default(std::string const & k, const char * x) const { return (defined(k) ? boost::python::extract<std::string>(_d[k]) : std::string(x)); }

  //-----------------------------------

  private:

  class cval {
   protected : 
    boost::python::dict dic;
    std::string key;
   public:
    cval (boost::python::object dic_, std::string key_): dic(),key(key_)  { dic = boost::python::extract<boost::python::dict>(dic_);}
    operator boost::python::object() const { return dic[key]; }
#define CAST_OPERATOR(T) operator T () const{ return (boost::python::extract<T>(dic[key])); }
     CAST_OPERATOR(short)
     CAST_OPERATOR(unsigned short)
     CAST_OPERATOR(int)
     CAST_OPERATOR(unsigned int)
     CAST_OPERATOR(long)
     CAST_OPERATOR(unsigned long)
     CAST_OPERATOR(long long)
     CAST_OPERATOR(unsigned long long)
     CAST_OPERATOR(float)
     CAST_OPERATOR(double)
     CAST_OPERATOR(long double)
     CAST_OPERATOR(bool)
     CAST_OPERATOR(std::string)
#undef CAST_OPERATOR
  };

  //-----------------------------------

  class val : public cval  { 
   public : 
    val (boost::python::object dic_, std::string key_): cval(dic_,key_) {}
    template <typename T> val & operator=(T const & rhs) { this->dic[this->key] = rhs; return *this; }
  };

  //-----------------------------------
  
  bool defined(std::string const & k) const { return _d.has_key(k); }

  void defined_or_throw(std::string const & k) const { 
   if (defined(k)) return;
   std::stringstream fs; 
   fs<<"class improved_python_dict : the key "<<k<<" is not present in the dictionnary";
   throw fs.str();
  }

 };

} } 

#endif

