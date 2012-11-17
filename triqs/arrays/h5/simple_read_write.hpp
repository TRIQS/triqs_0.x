/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011 by O. Parcollet
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
#include "./group_or_file.hpp"
#include <triqs/arrays/vector.hpp>
#include <string.h>

namespace triqs {
 namespace arrays {
  namespace h5 {

   /********************* IO for basic types and complex ... *************************************************/

   template <typename S> ENABLE_IF(boost::is_arithmetic<S>)
    h5_write (h5::group_or_file f, std::string const & name, S const & A) {
     try {
      f.unlink_if_exists(name);
      DataSet ds;
      ds = f->createDataSet( name.c_str(), data_type_file(S()), H5::DataSpace() );
      ds.write( (void *)(&A), data_type_mem_scalar(A), H5::DataSpace() );
     }
     TRIQS_ARRAYS_H5_CATCH_EXCEPTION;
    }

   template <typename S> ENABLE_IF(boost::is_arithmetic<S>)
    h5_read (h5::group_or_file f, std::string const & name,  S & A) {
     if (!f.exists(name))  TRIQS_RUNTIME_ERROR << "no such dataset : "<<name <<" in file ";
     try {
      DataSet ds = f->openDataSet( name.c_str() );
      DataSpace dataspace = ds.getSpace();
      int rank = dataspace.getSimpleExtentNdims();
      if (rank != 0) TRIQS_RUNTIME_ERROR << "triqs::array::h5::read. Rank mismatch : expecting a scalar (rank =0)"
       <<" while the array stored in the hdf5 file has rank = "<<rank;
      ds.read( (void *)(&A), data_type_mem_scalar(A), DataSpace() , DataSpace() );
     }
     TRIQS_ARRAYS_H5_CATCH_EXCEPTION;
    }

   /********************* IO for string  *************************************************/

   /**
    * \brief Write a string  into an hdf5 file
    * \param f The h5 file or group of type H5::H5File or H5::Group
    * \param name The name of the hdf5 array in the file/group where the stack will be stored
    * \param value The string
    * \exception The HDF5 exceptions will be caught and rethrown as TRIQS_RUNTIME_ERROR (with a full stackstrace, cf triqs doc).
    */
   inline void h5_write (h5::group_or_file f, std::string const & name, std::string const & value) {
    try {
     DataSet ds;
     StrType strdatatype(PredType::C_S1, value.size());
     ds = f->createDataSet( name.c_str(), strdatatype, DataSpace() );
     ds.write((void*)(value.c_str()), strdatatype );
    }
    TRIQS_ARRAYS_H5_CATCH_EXCEPTION;
   }

   /**
    * \brief Read a string from an hdf5 file
    * \param f The h5 file or group of type H5::H5File or H5::Group
    * \param name The name of the hdf5 array in the file/group where the stack will be stored
    * \param value The string to fill
    * \exception The HDF5 exceptions will be caught and rethrown as TRIQS_RUNTIME_ERROR (with a full stackstrace, cf triqs doc).
    */
   inline void h5_read (h5::group_or_file f, std::string const & name, std::string & value) {
    if (!f.exists(name))  TRIQS_RUNTIME_ERROR << "no such dataset : "<<name <<" in file ";
    try {
     DataSet ds = f->openDataSet( name.c_str() );
     DataSpace dataspace = ds.getSpace();
     int rank = dataspace.getSimpleExtentNdims();
     if (rank != 0) TRIQS_RUNTIME_ERROR << "Reading a string and got rank !=0";
     size_t size = ds.getStorageSize();
     StrType strdatatype(PredType::C_S1, size);
     std::vector<char> buf(size+1, 0x00);
     ds.read( (void *)(&buf[0]), strdatatype, DataSpace(), DataSpace() );
     value = ""; value.append( &(buf.front()) );
    }
    TRIQS_ARRAYS_H5_CATCH_EXCEPTION;
   }

   /*********************************** WRITE array ****************************************************************/
   /**
    * \brief Write an array or a view into an hdf5 file
    * \tparam ArrayType The type of the array/matrix/vector, etc..
    * \param f The h5 file or group of type H5::H5File or H5::Group
    * \param name The name of the hdf5 array in the file/group where the stack will be stored
    * \param A The array to be stored
    * \param C_reorder bool If true [default] the data will be stored in C order in the hdf5, hence making a temporary
    *        cache of the data to reorder them in memory.
    *        If false, the array is stored as it [if you know what you are doing]
    * \exception The HDF5 exceptions will be caught and rethrown as TRIQS_RUNTIME_ERROR (with a full stackstrace, cf triqs doc).
    */
   template <typename ArrayType>
    void write_array (group_or_file f, std::string const & name, ArrayType const & A, bool C_reorder = true) {
     try {
      if (h5::exists(f, name)) f->unlink( name.c_str());  // put some option here ?
      DataSet ds;
      if (C_reorder) {
       BOOST_AUTO(C, make_const_cache(A,Option::C()));
       ds = f->createDataSet( name.c_str(), data_type_file(typename ArrayType::value_type()), data_space(C.view()) );
       ds.write( data(C.view()), data_type_mem(A), data_space(C.view()) );
      }
      else {
       ds = f->createDataSet( name.c_str(), data_type_file(typename ArrayType::value_type()), data_space(A) );
       ds.write( data(A), data_type_mem(A), data_space(A) );
      }
      // if complex, to be python compatible, we add the __complex__ attribute
      if (boost::is_complex<typename ArrayType::value_type>::value)  write_attribute(ds,"__complex__","1");
     }
     TRIQS_ARRAYS_H5_CATCH_EXCEPTION;
    }

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
    void read_array (group_or_file f, std::string const & name,  ArrayType & A, bool C_reorder = true) {
     typedef typename ArrayType::value_type V;
     if (!h5::exists(f, name))  TRIQS_RUNTIME_ERROR << "no such dataset : "<<name <<" in file ";
     try {
      DataSet ds = f->openDataSet( name.c_str() );
      DataSpace dataspace = ds.getSpace();
      static const unsigned int Rank =  ArrayType::rank + (boost::is_complex<typename ArrayType::value_type>::value ? 1 : 0);
      int rank = dataspace.getSimpleExtentNdims();
      if (rank != Rank) TRIQS_RUNTIME_ERROR << "triqs::array::h5::read. Rank mismatch : the array has rank = "
       <<Rank<<" while the array stored in the hdf5 file has rank = "<<rank;
      mini_vector<hsize_t,Rank> dims_out;
      //int ndims = dataspace.getSimpleExtentDims( &dims_out[0], NULL);
      dataspace.getSimpleExtentDims( &dims_out[0], NULL);
      mini_vector<size_t,ArrayType::rank > d2; for (size_t u=0; u<ArrayType::rank ; ++u) d2[u] = dims_out[u];
      resize_or_check(A, d2 );
      if (C_reorder) {
       BOOST_AUTO(C,  make_cache(A, Option::C() ));
       ds.read( data(C.view()), data_type_mem(C.view()), data_space(C.view()) , dataspace );
      }
      else { ds.read( data(A), data_type_mem(A), data_space(A) , dataspace ); }
     }
     TRIQS_ARRAYS_H5_CATCH_EXCEPTION;
    }

   /*********************************** 1d array and std::vector of string *********************************/

   template<typename ArrayVectorOfStringType>
    void write_1darray_vector_of_string_impl (group_or_file f, std::string const & name, ArrayVectorOfStringType const & V) {
     size_t s=0; for (auto & x : V) s = std::max(s,x.size());
     try {
      if (h5::exists(f, name)) f->unlink( name.c_str());  // put some option here ?
      DataSet ds;
      StrType strdatatype(PredType::C_S1, s);
      const size_t n = V.size();
      std::vector<char> buf(n*s, 0x00);
      size_t i=0; for (auto & x : V) {strcpy( &buf[i*s], x.c_str()); ++i;}

      mini_vector<hsize_t,1> L; L[0]=V.size();
      mini_vector<hsize_t,1> S; S[0]=1;
      auto d_space = dataspace_from_LS<1,false > (L,L,S);

      ds = f->createDataSet( name.c_str(), strdatatype, d_space );
      ds.write( (void *)(&buf[0]),strdatatype, d_space );
     }
     TRIQS_ARRAYS_H5_CATCH_EXCEPTION;
    }

   template<typename ArrayVectorOfStringType>
    void read_1darray_vector_of_string_impl (group_or_file f, std::string const & name, ArrayVectorOfStringType & V) {
     if (!h5::exists(f, name))  TRIQS_RUNTIME_ERROR << "no such dataset : "<<name <<" in file ";
     try {
      DataSet ds = f->openDataSet( name.c_str() );
      DataSpace dataspace = ds.getSpace();
      mini_vector<hsize_t,1> dims_out;
      int ndims = dataspace.getSimpleExtentDims( &dims_out[0], NULL);
      if (ndims !=1) TRIQS_RUNTIME_ERROR << "triqs::h5 : Trying to read 1d array/vector . Rank mismatch : the array stored in the hdf5 file has rank = "<<ndims;

      size_t Len = dims_out[0];
      V.resize(Len);
      size_t size = ds.getStorageSize();
      StrType strdatatype(PredType::C_S1, size);

      std::vector<char> buf(Len*(size+1), 0x00);

      mini_vector<hsize_t,1> L; L[0]=V.size();
      mini_vector<hsize_t,1> S; S[0]=1;
      auto d_space = dataspace_from_LS<1,false > (L,L,S);

      ds.read( (void *)(&buf[0]),strdatatype, d_space );
      size_t i=0; for (auto & x : V) { x = ""; x.append(&buf[i*(size)]); ++i;}
     }
     TRIQS_ARRAYS_H5_CATCH_EXCEPTION;
    }
  }// end of h5 namespace



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
   h5_read (h5::group_or_file fg, std::string const & name,  ArrayType & A) { h5::read_array(fg,name, A);}

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
   h5_write (h5::group_or_file fg, std::string const & name,  ArrayType const & A) { h5::write_array(fg,name, A);}


  inline void h5_write (h5::group_or_file f, std::string const & name, vector_view<std::string> const & V) {
   h5::write_1darray_vector_of_string_impl(f,name,V);
  }

  inline void h5_read (h5::group_or_file f, std::string const & name, arrays::vector<std::string> & V) {
   h5::read_1darray_vector_of_string_impl(f,name,V);
  }

  // I can not use the generic code, just because the resize of the array take a shape,  not a size_t as std::vector and vector
  // Ok, speed is no issue here...
  inline void h5_read (h5::group_or_file f, std::string const & name, arrays::array<std::string,1> & V) {
   arrays::vector<std::string> res; h5_read(f,name,res); V = res;
  }

 }
}

// for ADL, need to put this in std:: namespace ...
namespace std {
 inline void h5_write (triqs::arrays::h5::group_or_file f, std::string const & name, std::vector<std::string> const & V) {
  triqs::arrays::h5::write_1darray_vector_of_string_impl(f,name,V);
 }

 inline void h5_read (triqs::arrays::h5::group_or_file f, std::string const & name, std::vector<std::string> & V) {
  triqs::arrays::h5::read_1darray_vector_of_string_impl(f,name,V);
 }
}
#endif

