/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011-2013 by O. Parcollet
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
#ifndef TRIQS_ARRAYS_H5_LOWLEVEL_H
#define TRIQS_ARRAYS_H5_LOWLEVEL_H
#include <triqs/arrays/array.hpp>
#include <triqs/arrays/vector.hpp>
#include <triqs/h5.hpp>
#include "../cache.hpp"

namespace triqs { namespace arrays { 
 namespace h5_impl {

  template <typename T, int R> const void * get_data_cptr ( array_view<T,R> const & A) { return h5::get_data_ptr(&(A.storage()[0]));}
  template <typename T, int R> const void * get_data_cptr ( array<T,R> const & A) { return h5::get_data_ptr(&(A.storage()[0]));}
  template <typename T, int R> void * get_data_ptr ( array_view<T,R> & A) { return h5::get_data_ptr(&(A.storage()[0]));}
  template <typename T, int R> void * get_data_ptr ( array<T,R> & A) { return h5::get_data_ptr(&(A.storage()[0]));}
 
  // the dataspace corresponding to the array. Contiguous data only...
  template <typename ArrayType >
   H5::DataSpace data_space ( ArrayType const & A) { 
    static const unsigned int R = ArrayType::rank;
    mini_vector<hsize_t,R> S; 
    mini_vector<std::ptrdiff_t,R> const & S1 ( A.indexmap().strides() );
    for (size_t u=0; u<R ; ++u) { 
     if (S1[u]<=0) TRIQS_RUNTIME_ERROR<<" negative strides not permitted in h5"; 
     S[u] =1; 
    }   
    if (!A.indexmap().is_contiguous())  TRIQS_RUNTIME_ERROR<<" h5 : internal error : array not contiguous";
    static const bool is_complex =  boost::is_complex<typename ArrayType::value_type>::value;
    return h5::dataspace_from_LS<R,is_complex > ( A.indexmap().domain().lengths(),A.indexmap().domain().lengths(), S);
   }

  /********************   resize or check the size ****************************************************/

  template <typename A> ENABLE_IF(is_amv_value_class<A>) 
   resize_or_check ( A & a, mini_vector<size_t,A::rank> const & dimsf ) { a.resize( indexmaps::cuboid::domain_t<A::rank>( dimsf)); }

  template <typename A> ENABLE_IF(is_amv_view_class<A>) 
   resize_or_check ( A const & a, mini_vector<size_t,A::rank> const & dimsf ) { 
    if (a.indexmap().domain().lengths() != dimsf) TRIQS_RUNTIME_ERROR<<"Dimension error : the view can not be resized : " 
     << "\n in file  : "<< dimsf.to_string() 
      << "\n in view  : "<<a.indexmap().domain().lengths().to_string() ;
   }

  /*********************************** WRITE array ****************************************************************/
  /**
   * \brief Write an array or a view into an hdf5 file
   * \tparam 
   * \param f The h5 file or group of type H5::H5File or H5::Group
   * \param name The name of the hdf5 array in the file/group where the stack will be stored
   * \param A The array to be stored
   * \param C_reorder bool If true [default] the data will be stored in C order in the hdf5, hence making a temporary
   *        cache of the data to reorder them in memory.
   *        If false, the array is stored as it [if you know what you are doing]
   * \exception The HDF5 exceptions will be caught and rethrown as TRIQS_RUNTIME_ERROR (with a full stackstrace, cf triqs doc).
   */
  template <typename T, int R>
   void write_array (h5::group g, std::string const & name, array_view <T,R> const & A, bool C_reorder = true) {
    if (C_reorder) { write_array(g,name, make_const_cache(A).view(),false);}
    try {
     H5::DataSet ds = g.create_dataset(name, h5::data_type_file<T>(), data_space(A) );
     ds.write( get_data_cptr(A), h5::data_type_memory<T>(), data_space(A) );
     // if complex, to be python compatible, we add the __complex__ attribute
     if (boost::is_complex<T>::value)  h5::write_string_attribute(&ds,"__complex__","1");
    }
    TRIQS_ARRAYS_H5_CATCH_EXCEPTION;
   }
 
  template <typename T, int R>
   void write_array (h5::group g, std::string const & name, array <T,R> const & A, bool C_reorder = true) { write_array(g,name,A(),C_reorder);}
 
  /*********************************** READ array ****************************************************************/
  /**
   * \brief Read an array or a view from an hdf5 file
   * \tparam ArrayType The type of the array/matrix/vector, etc..
   * \param f The h5 file or group of type H5::H5File or H5::Group
   * \param name The name of the hdf5 array in the file/group where the stack will be stored
   * \param A The array to be stored
   * \param C_reorder bool If true [default] the data will be stored in C order in the hdf5, hence making a temporary
   *        cache of the data to reorder them in memory. If false, the array is stored as it [if you know what you are doing]
   * \exception The HDF5 exceptions will be caught and rethrown as TRIQS_RUNTIME_ERROR (with a full stackstrace, cf triqs doc).
   */
  template <typename ArrayType>
   void read_array (h5::group g, std::string const & name,  ArrayType & A, bool C_reorder = true) {
    typedef typename ArrayType::value_type V;
    try {
     H5::DataSet ds = g.open_dataset(name);
     H5::DataSpace dataspace = ds.getSpace();
     static const unsigned int Rank =  ArrayType::rank + (boost::is_complex<typename ArrayType::value_type>::value ? 1 : 0);
     int rank = dataspace.getSimpleExtentNdims();
     if (rank != Rank) TRIQS_RUNTIME_ERROR << "triqs::array::h5::read. Rank mismatch : the array has rank = "
      <<Rank<<" while the array stored in the hdf5 file has rank = "<<rank;
     mini_vector<hsize_t,Rank> dims_out;
     //int ndims = dataspace.getSimpleExtentDims( &dims_out[0], NULL);
     dataspace.getSimpleExtentDims( &dims_out[0], NULL);
     mini_vector<size_t,ArrayType::rank > d2; for (size_t u=0; u<ArrayType::rank ; ++u) d2[u] = dims_out[u];
     resize_or_check(A, d2 );
     if (C_reorder) { read_array(g,name, make_cache(A).view(),false);}
     ds.read( get_data_ptr(A), h5::data_type_memory<typename ArrayType::value_type>(), data_space(A) , dataspace );
    }
    TRIQS_ARRAYS_H5_CATCH_EXCEPTION;
   }

 }// namespace

 template<typename ArrayType> struct is_amv_value_or_view_class_no_string :
  boost::mpl::and_<is_amv_value_or_view_class<ArrayType>, boost::mpl::not_<boost::is_base_of<std::string, typename ArrayType::value_type> > > {};

 /**
  * \brief Read an array or a view from an hdf5 file
  * \tparam ArrayType The type of the array/matrix/vector, etc..
  * \param fg The h5 file or group of type H5::H5File or H5::Group
  * \param name The name of the hdf5 array in the file/group where the stack will be stored
  * \param A The array to be stored
  * \exception The HDF5 exceptions will be caught and rethrown as TRIQS_RUNTIME_ERROR (with a full stackstrace, cf triqs doc).
  */
 template <typename ArrayType>
  ENABLE_IF(is_amv_value_or_view_class_no_string<ArrayType>)
  h5_read (h5::group fg, std::string const & name,  ArrayType & A) { h5_impl::read_array(fg,name, A);}

 /**
  * \brief Write an array or a view into an hdf5 file
  * \tparam ArrayType The type of the array/matrix/vector, etc..
  * \param fg The h5 file or group of type H5::H5File or H5::Group
  * \param name The name of the hdf5 array in the file/group where the stack will be stored
  * \param A The array to be stored
  * \exception The HDF5 exceptions will be caught and rethrown as TRIQS_RUNTIME_ERROR (with a full stackstrace, cf triqs doc).
  */
 template <typename ArrayType>
  ENABLE_IF(is_amv_value_or_view_class_no_string<ArrayType>)
  h5_write (h5::group fg, std::string const & name,  ArrayType const & A) { h5_impl::write_array(fg,name, A);}


 inline void h5_write (h5::group f, std::string const & name, vector_view<std::string> const & V) {
  h5::detail::write_1darray_vector_of_string_impl(f,name,V);
 }

 inline void h5_read (h5::group f, std::string const & name, arrays::vector<std::string> & V) {
  h5::detail::read_1darray_vector_of_string_impl(f,name,V);
 }

 // I can not use the generic code, just because the resize of the array take a shape,  not a size_t as std::vector and vector
 // Ok, speed is no issue here...
 inline void h5_read (h5::group f, std::string const & name, arrays::array<std::string,1> & V) {
  arrays::vector<std::string> res; h5_read(f,name,res); V = res;
 }

}}
#endif

