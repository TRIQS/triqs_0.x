/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012 by M. Ferrero, O. Parcollet
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
#ifndef TRIQS_GF_TOOLS_H
#define TRIQS_GF_TOOLS_H
#include <triqs/utility/first_include.hpp>
#include <utility>
#include <boost/iterator/iterator_facade.hpp>
#include <triqs/utility/cint.hpp>
#include <triqs/clef/core.hpp>
#include <triqs/arrays/array.hpp>
#include <triqs/arrays/matrix.hpp>
#include <triqs/utility/proto/tools.hpp>
#include <triqs/arrays/h5/simple_read_write.hpp>
#include <triqs/arrays/proto/matrix_algebra.hpp>
#include <triqs/arrays/proto/array_algebra.hpp>
#include "triqs/utility/complex_ops.hpp"
#include <triqs/utility/view_tools.hpp>

namespace triqs { namespace gf {

 namespace tqa= triqs::arrays;
 namespace mpl=boost::mpl;
 namespace proto=boost::proto;
 namespace bf = boost::fusion; 
 namespace tup = triqs::utility::proto;

 //------------------------------------------------------

 typedef std::complex<double> dcomplex; 

 enum statistic_enum {Boson,Fermion};

 struct freq_infty{}; // the point at infinity

 inline  std::vector<std::string> split(const std::string &s, char delim){
  std::vector<std::string> elems;
  std::stringstream ss(s);
  std::string item;
  while(std::getline(ss, item, delim)) { elems.push_back(item); }
  return elems;
 }

 //------------------------------------------------------

 class indices_2_t { 
  std::vector<std::vector<std::string>> data;
  public : 
  indices_2_t() {}
  indices_2_t(std::vector<std::vector<std::string>> const & d) : data(d) {} 
  
  indices_2_t(int n1, int n2) { 
   std::vector<std::string> v1,v2;
   for (int i =0; i<n1; ++i)  { std::stringstream fs; fs<<i; v1.push_back(fs.str()); }
   for (int i =0; i<n2; ++i)  { std::stringstream fs; fs<<i; v2.push_back(fs.str()); }
   data.push_back(v1); data.push_back(v2);
  }
 
  std::vector<std::string> const & operator[](int i) const { return data[i];} 
  std::vector<std::vector<std::string>> const & operator()() const { return data;} 

  bool same() const { return data[0]==data[1];} 

  /// Write into HDF5
  friend void h5_write (tqa::h5::group_or_file fg, std::string key, indices_2_t const & g) {
   h5_write(fg,key,g.data[0]);
  }

  /// Read from HDF5
  friend void h5_read  (tqa::h5::group_or_file fg, std::string key, indices_2_t & g){
   std::vector<std::string> V;
   h5_read(fg,key,V);
   g.data.push_back(V); g.data.push_back(V);
  }

  //  BOOST Serialization
  friend class boost::serialization::access;
  template<class Archive>
   void serialize(Archive & ar, const unsigned int version) {
    ar & boost::serialization::make_nvp("data",data);
   }

 };

 //------------------------------------------------------

 struct nothing {
  template<typename... Args> nothing(Args...) {} // takes anything, do nothing..
  nothing() {}
  typedef void has_view_type_tag;     // Idiom : ValueView  
  typedef nothing view_type;
  typedef nothing non_view_type;
  void rebind (nothing){}
  void operator=(nothing) {}
  friend void h5_write (tqa::h5::group_or_file fg, std::string subgroup_name, nothing ) {}
  friend void h5_read  (tqa::h5::group_or_file fg, std::string subgroup_name, nothing ) {}
  //  BOOST Serialization
  friend class boost::serialization::access;
  template<class Archive>
   void serialize(Archive & ar, const unsigned int version) {
   }
  friend nothing operator +( nothing, nothing) { return nothing();}
 }; 

 //------------------------------------------------------

}}
#endif
