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
#include "./common.hpp"
namespace triqs { 
 namespace arrays { 
  namespace h5 { 

   /*********************************** WRITE array ****************************************************************/
   /** 
    * \brief Write an array or a view into an hdf5 file
    * \tparam ArrayType The type of the array/matrix/vector, etc..
    * \tparam FileGroupType The type of the file or of the group
    * \param file_or_group The h5 file or group, of type FileGroupType
    * \param name The name of the hdf5 array in the file/group where the stack will be stored
    * \param A The array to be stored
    * \param C_reorder bool If true [default] the data will be stored in C order in the hdf5, hence making a temporary
    *        cache of the data to reorder them in memory. 
    *        If false, the array is stored as it [if you know what you are doing]
    * \exception The HDF5 exceptions will be caught and rethrown as TRIQS_RUNTIME_ERROR (with a full stackstrace, cf triqs doc). 
    */
   template <typename ArrayType, typename FileGroupType >
    void write (FileGroupType file_or_group, std::string const & name, ArrayType const & A, bool C_reorder = true) { 
     try { 
      if (h5::exists(file_or_group, name)) file_or_group.unlink( name.c_str());  // put some option here ?
      DataSet ds;
      if (C_reorder) { 
       BOOST_AUTO(C, make_const_cache(A,Option::C()));
       //typename result_of::cache<false,Tag::C, ArrayType >::type C(A);
       ds = file_or_group.createDataSet( name.c_str(), data_type_file(typename ArrayType::value_type()), data_space(C.view()) );
       ds.write( data(C.view()), data_type_mem(A), data_space(C.view()) ); 
      }
      else { 
       ds = file_or_group.createDataSet( name.c_str(), data_type_file(typename ArrayType::value_type()), data_space(A) );
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
    * \tparam FileGroupType The type of the file or of the group
    * \param file_or_group The h5 file or group, of type FileGroupType
    * \param name The name of the hdf5 array in the file/group where the stack will be stored
    * \param A The array to be stored
    * \param C_reorder bool If true [default] the data will be stored in C order in the hdf5, hence making a temporary
    *        cache of the data to reorder them in memory. If false, the array is stored as it [if you know what you are doing]
    * \exception The HDF5 exceptions will be caught and rethrown as TRIQS_RUNTIME_ERROR (with a full stackstrace, cf triqs doc). 
    */
   template <typename ArrayType, typename FileGroupType >
    void read (FileGroupType file_or_group, std::string const & name,  ArrayType & A, bool C_reorder = true) { 
     typedef typename ArrayType::value_type V;
     if (!h5::exists(file_or_group, name))  TRIQS_RUNTIME_ERROR << "no such dataset : "<<name <<" in file ";
     try { 
      DataSet ds = file_or_group.openDataSet( name.c_str() );
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
       //typename result_of::cache<true,Tag::C, ArrayType >::type C(A);
       ds.read( data(C.view()), data_type_mem(C.view()), data_space(C.view()) , dataspace );
      }
      else { ds.read( data(A), data_type_mem(A), data_space(A) , dataspace ); } 
     } 
     TRIQS_ARRAYS_H5_CATCH_EXCEPTION;
    }
  }

  template <typename ArrayType, typename FileGroupType >
   typename boost::enable_if<arrays::is_amv_value_or_view_class<ArrayType> >::type 
   h5_read (FileGroupType file_or_group, std::string const & name,  ArrayType & A) { h5::read(file_or_group,name, A);} 

  template <typename ArrayType, typename FileGroupType >
   typename boost::enable_if<arrays::is_amv_value_or_view_class<ArrayType> >::type 
   h5_write (FileGroupType file_or_group, std::string const & name,  ArrayType const & A) { h5::write(file_or_group,name, A);} 

 }}
#endif

