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
/// triqs  
namespace triqs {
 /// arrays
 namespace arrays { 
  /// h5
  namespace h5 { 

 /*********************************** WRITE array ****************************************************************/

 /** 
  * Write an array or a view in an hdf5 file
  */
 template <typename ArrayType, typename FileGroupType >
  void write (FileGroupType file_or_group, std::string const & name, ArrayType const & A, bool C_reorder = true) { 
   try { 
    if (h5::exists(file_or_group, name)) file_or_group.unlink( name.c_str());  // put some option here ?
    DataSet ds;
    if (C_reorder) { 
     typename result_of::cache<false,Tag::C, ArrayType >::type C(A);
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
  * Read an array or a view from an hdf5 file
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
    if (rank != Rank) TRIQS_RUNTIME_ERROR << "dimension error : expected "<<Rank<<" got "<<rank;
    mini_vector<hsize_t,Rank> dims_out;
    int ndims = dataspace.getSimpleExtentDims( &dims_out[0], NULL);
    mini_vector<size_t,ArrayType::rank > d2; for (size_t u=0; u<ArrayType::rank ; ++u) d2[u] = dims_out[u];
    resize_or_check(A, d2 ); 
    if (C_reorder) { 
     typename result_of::cache<true,Tag::C, ArrayType >::type C(A);
     ds.read( data(C.view()), data_type_mem(C.view()), data_space(C.view()) , dataspace );
    }
    else { ds.read( data(A), data_type_mem(A), data_space(A) , dataspace ); } 
   } 
   TRIQS_ARRAYS_H5_CATCH_EXCEPTION;
  }
}}}
#endif

