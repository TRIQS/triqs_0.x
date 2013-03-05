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
#ifndef TRIQS_GF_DATA_PROXIES_H
#define TRIQS_GF_DATA_PROXIES_H
#include <triqs/utility/first_include.hpp>
#include <utility>
#include <triqs/arrays.hpp>

namespace triqs { namespace gf {

 template<typename T, int R> struct data_proxy_array;

 template<typename T> struct data_proxy_array<T,3> { 
  /// The storage
  typedef arrays::array<T,3>             storage_t;
  typedef typename storage_t::view_type  storage_view_t;

  /// The data access 
  arrays::view_proxy<storage_t,2>       operator()(storage_t       & data, size_t i)       { return arrays::view_proxy<storage_t,2>(data,i); } 
  arrays::const_view_proxy<storage_t,2> operator()(storage_t const & data, size_t i) const { return arrays::const_view_proxy<storage_t,2>(data,i); } 
  arrays::view_proxy<storage_view_t,2>       operator()(storage_view_t       & data, size_t i)       { return arrays::view_proxy<storage_view_t,2>(data,i); } 
  arrays::const_view_proxy<storage_view_t,2> operator()(storage_view_t const & data, size_t i) const { return arrays::const_view_proxy<storage_view_t,2>(data,i); } 
 };

 template<typename T> struct data_proxy_array<T,1>{ 
  /// The storage
  typedef arrays::array<T,1>             storage_t;
  typedef typename storage_t::view_type  storage_view_t;

  /// The data access 
  auto operator()(storage_t       & data,size_t i)       -> decltype(data(i)) { return data(i);}
  auto operator()(storage_t const & data,size_t i) const -> decltype(data(i)) { return data(i);}
  auto operator()(storage_view_t       & data,size_t i)       -> decltype(data(i)) { return data(i);}
  auto operator()(storage_view_t const & data,size_t i) const -> decltype(data(i)) { return data(i);}
 };

 template<typename T> struct data_proxy_vector { 
  typedef typename T::view_type Tv;

  /// The storage
  typedef std::vector<T> storage_t;
  typedef std::vector<Tv> storage_view_t;

  /// The data access 
  T       &  operator()(std::vector<T> &       data, size_t i)  { return data[i];}
  T const &  operator()(std::vector<T> const & data, size_t i)  { return data[i];}
  Tv       &  operator()(std::vector<Tv> &       data, size_t i)  { return data[i];}
  Tv const &  operator()(std::vector<Tv> const & data, size_t i)  { return data[i];}
 };

}}
#endif

