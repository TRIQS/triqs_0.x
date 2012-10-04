
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

#ifndef TRIQS_GF_C_H 
#define TRIQS_GF_C_H

#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <list>
#include <vector>
#include <triqs/utility/report_stream.hpp>
#include <triqs/utility/exceptions.hpp>
#include <algorithm> 
#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp> 

/** 
  This is a little wrapper to see the python class GF in C++
  GFBLOC must be the type of the GF blocks
  */
template<class GFBLOC>
class GF_C 
{
 public:

  typedef std::pair<std::string, boost::shared_ptr<GFBLOC> > ElementType ;
  typedef std::vector<ElementType>  vector_type;

 protected:

  vector_type data;
  static bool cmp(ElementType const &a, ElementType const & b) { return (a.second->N1 < b.second->N1);}
  inline size_t compute_block_size_max() { return std::max_element(data.begin(), data.end(), &GF_C<GFBLOC>::cmp)->second->N1;}
  size_t _nblocks, _block_size_max;
  double _beta;

 public: 

  int NBlocks() const { return _nblocks;}
  int max_block_size () const { return _block_size_max;}
  double Beta () const { return _beta;}

  GF_C(vector_type const & _v): data(_v) { 
   for (typename vector_type::const_iterator it=data.begin(); it !=data.end(); ++it)  assert(it->second);
   _block_size_max =  compute_block_size_max();
   _nblocks  = data.size();
   _beta = data[0].second->Beta;
  }

  // returns a deep copy of the object
  GF_C deep_copy() const { 
   vector_type v; 
   for (typename vector_type::const_iterator it=data.begin(); it !=data.end(); ++it) 
    v.push_back( std::make_pair(it->first, boost::make_shared<GFBLOC> (*(it->second))));
   return GF_C(v);
  }

 private: //not implemented
  GF_C<GFBLOC> & operator = (const GF_C<GFBLOC> & Gin);

 public :

  GFBLOC & operator[](int bl) {  assert ((bl< int(_nblocks))); assert (data[bl].second); return( *(data[bl].second)) ;}
   const GFBLOC & operator[] (int bl) const { assert ((bl<int(_nblocks))); return( *(data[bl].second)) ;}

   int block_size (int bl) const {return (*this)[bl].N1;}  

 /*  inline void operator +=  (const GF_C<GFBLOC> & Gin) 
   { if (!same_structure(Gin)) TRIQS_RUNTIME_ERROR<<"Gf -= : structure is not the same"; for (int bl =0; bl<_nblocks; bl++) (*this)[bl] += Gin[bl];}

   inline void operator -=  (const GF_C<GFBLOC> & Gin) 
   { if (!same_structure(Gin)) TRIQS_RUNTIME_ERROR<<"Gf -= : structure is not the same"; for (int bl =0; bl<_nblocks; bl++) (*this)[bl] -= Gin[bl];}

   inline void operator *= (const GF_C<GFBLOC> & Gin) { for (int bl =0; bl<_nblocks; bl++) (*this)[bl] *= Gin[bl];}
   inline void operator *= ( double alpha) { for (int bl =0; bl<_nblocks; bl++) (*this)[bl] *= alpha; }
   template<typename T>
    inline void operator /=  (const T & alpha)  { for (int bl =0; bl<_nblocks; bl++) (*this)[bl] /= alpha;}
*/
   //inline void save(std::string file,  bool accumulate) { boost::python::call_method<boost::python::object,std::string,bool>(this->ptr(),"save",file,accumulate);}
   //inline void load(std::string file) { boost::python::call_method<boost::python::object,std::string>(this->ptr(),"load",file);}

   void zero() {for (int bl =0; bl<_nblocks; bl++) (*this)[bl].zero();}

 protected:  
   //   Check that the structure of the 2 Greens function is the same
   inline bool same_structure (const GF_C<GFBLOC> & Gin) const {
    bool res = (Gin._nblocks == this->_nblocks);
    if (!res) return res;
    for (int i=0; i<this->_nblocks; ++i)  res = res && (Gin[i].N1 == (*this)[i].N1);
    return res;
   }

};

namespace GF_C_details {
 using namespace boost;

 inline PyObject* to_pyo( python::object const & obj) { 
  python::object Res =  obj;// keep it until the end of the routine
  PyObject * res = Res.ptr();
  Py_INCREF(res); return res;
 }  

 template<class GFBLOC>
  struct GF_C_to_python  { 
   static PyObject* convert (GF_C<GFBLOC>  const & ) { 
    assert(0);
    python::object obj; // TO BE WRITTEN 
    return to_pyo(obj); 
   }  
  };


 // simply use the python_tools converter !
 template<class GFBLOC> 
  struct GF_C_from_python_str {

   typedef GF_C<GFBLOC> Type;
   typedef typename Type::vector_type GF_stack_type;

   GF_C_from_python_str(){ python::converter::registry::push_back(&convertible, &construct,python::type_id< Type >()); }

   static GF_stack_type make_GF_Stack(boost::python::object OB) {
    GF_stack_type v;
    boost::python::stl_input_iterator<boost::python::object> begin_list(OB), end_list;
    std::list<boost::python::object>  l(begin_list,end_list);
    for (std::list<boost::python::object>::const_iterator p = l.begin(); p !=l.end(); ++p) { 
     std::string name = boost::python::extract<std::string> ((*p)[0]);
     boost::shared_ptr<GFBLOC> g = boost::python::extract<boost::shared_ptr<GFBLOC> > ((*p)[1]);
     if (g->N1!=g->N2) TRIQS_RUNTIME_ERROR<<"make_GF_Stack : non square matrix !";
     v.push_back(std::make_pair(name,g));
    }
    return v;
   }

   static void* convertible(PyObject* obj_ptr) {
    python::object obj ( python::borrowed (obj_ptr));
    bool ok = true;
    try  {
     boost::python::stl_input_iterator<boost::python::object> begin_list(obj), end_list;
     std::list<boost::python::object>  l(begin_list,end_list);
     for (std::list<boost::python::object>::const_iterator p = l.begin(); p !=l.end(); ++p) { 
      boost::python::extract<GFBLOC> ch((*p)[1]); // check that p is a tuple of 2-tuple...
      ok = ok && ch.check();
     } 
    }
    catch(...) { ok=false;}
    return ( ok ? obj_ptr : 0 );
   } 

   static void construct( PyObject* obj_ptr, python::converter::rvalue_from_python_stage1_data* data) {
    typedef python::converter::rvalue_from_python_storage< Type > storage_t;
    storage_t* the_storage = reinterpret_cast<storage_t*>( data );
    void* memory_chunk = the_storage->storage.bytes;
    python::object obj ( python::borrowed (obj_ptr));
    GF_stack_type st (make_GF_Stack(obj));
    new (memory_chunk) Type(st);
    data->convertible = memory_chunk;
   }
  };

 template<class GFBLOC>  
  inline void register_converter() { 
   boost::python::to_python_converter< GF_C<GFBLOC>, GF_C_to_python<GFBLOC>    >();
   GF_C_from_python_str<GFBLOC> ();
  }  

} // GF_C_details

#endif
