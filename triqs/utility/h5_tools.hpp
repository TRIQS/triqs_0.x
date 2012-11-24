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
#ifndef TRIQS_ARRAYS_H5_TOOLS_H
#define TRIQS_ARRAYS_H5_TOOLS_H
#include <H5Cpp.h>
#include "./view_tools.hpp"
#include <string> 

namespace triqs {  namespace h5 {
 using namespace H5;

 /**
  *  \brief Opaque object to handle an H5 group or file
  */
 class group_or_file {
  boost::shared_ptr<H5::CommonFG> fg_ptr;
  hid_t id;
  public:
  /**
   * \brief Constructor
   *  \param g H5 group
   */
  group_or_file (H5::Group g)   { fg_ptr = boost::make_shared<H5::Group>(g); id = g.getId();}
  /**
   * \brief Constructor
   *  \param f H5 file
   */
  group_or_file (H5::H5File f) { fg_ptr = boost::make_shared<H5::H5File>(f); id = f.getId();}

  group_or_file() {}

  group_or_file(std::string filename, int flag) {
   H5::H5File f(filename, H5F_ACC_TRUNC);
   fg_ptr = boost::make_shared<H5::H5File>(f); id = f.getId();
  }

  group_or_file(hid_t id_, bool is_group) {
   if (is_group) {
    H5::Group g;
    g.setId(id_);
    fg_ptr = boost::make_shared<H5::Group>(g);
   }
   else {
    H5::H5File f;
    f.setId(id_);
    fg_ptr = boost::make_shared<H5::H5File>(f);
   }
   id = id_;
  }

  hid_t getId() const { return id;}
  ///
  const H5::CommonFG & operator* () const  { return *fg_ptr;}
  ///
  H5::CommonFG & operator* ()              { return *fg_ptr;}
  ///
  H5::CommonFG const * operator-> () const { return fg_ptr.get();}
  ///
  H5::CommonFG * operator-> ()             { return fg_ptr.get();}
  ///
  bool exists (std::string path) const { return (H5Lexists(getId(), path.c_str(), H5P_DEFAULT)); }
  ///
  void unlink_if_exists(std::string path) const { if (exists(path)) fg_ptr->unlink(path.c_str()); }
  /// Open a subgroup
  group_or_file open_group(std::string Name) const { return fg_ptr->openGroup(Name.c_str());}
  /// Create a subgroup (unlink an existing one).
  group_or_file create_group(std::string Name) const {
   unlink_if_exists(Name); // CHOICE : WHAT DO WE DECIDE AS A POLICY ???
   return fg_ptr->createGroup(Name.c_str());
  }
 };

 /********************* exists *************************************************/

 inline bool exists (group_or_file f, std::string name) { return f.exists(name); }

 /********************* extractor *************************************************/

 template<typename T> 
  struct h5_extractor { 
   T const & operator() (group_or_file fg, std::string subgroup_name) const { 
    typename non_view_type_if_exists_else_type<T>::type r;
    h5_read(fg,subgroup_name,r);
    return r;
   }
  };

 /********************* macro *************************************************/

#define TRIQS_ARRAYS_H5_CATCH_EXCEPTION \
 catch( triqs::runtime_error error)  { throw triqs::runtime_error() << error.what();}\
 catch( FileIException error ) { error.printError(); TRIQS_RUNTIME_ERROR<<"H5 File error"; }\
 catch( DataSetIException error ) { error.printError(); TRIQS_RUNTIME_ERROR<<"H5 DataSet error"; }\
 catch( DataSpaceIException error ) { error.printError();  TRIQS_RUNTIME_ERROR<<"H5 DataSpace error"; }\
 catch( DataTypeIException error ) { error.printError(); TRIQS_RUNTIME_ERROR<<"H5 DataType error";  }\
 catch(...) { TRIQS_RUNTIME_ERROR<<"H5 unknown error";} 

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

}}
#endif

