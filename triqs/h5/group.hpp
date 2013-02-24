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
#ifndef TRIQS_H5_GROUP_H
#define TRIQS_H5_GROUP_H
#include "./base.hpp"
namespace triqs { namespace h5 {
 /**
  *  \brief A local derivative of Group.
  *  Rational : use ADL for h5_read/h5_write, catch and rethrow exception, add some policy for opening/creating 
  */
 class group : public H5::Group {
  public:
   group() = default;
   group( group const & ) = default;
   group (H5::Group g) : H5::Group(g) {} 
   
   /// Takes the "/" group at the top of the file.
   group (H5::H5File f) : H5::Group(f.openGroup("/")) {} // can not fail, right ?

   group(std::string const & filename, int flag) { H5::H5File f(filename, H5F_ACC_TRUNC); *this = group(f);}

   group(hid_t id_, bool is_group) {
    if (is_group) { this->setId(id_); }
    else { H5::H5File f; f.setId(id_); *this = group(f); }
   }
   ///  Write the triqs tag of the group if it is an object.
   template<typename T>  void write_triqs_hdf5_data_scheme (T const & obj) {  
     write_string_attribute( this, "TRIQS_HDF5_data_scheme" , get_triqs_hdf5_data_scheme(obj).c_str() ) ;
    }   
   /// Read the triqs tag of the group if it is an object. Returns"" if attribute is not present
   std::string read_triqs_hdf5_data_scheme() const { return read_string_attribute(this, "TRIQS_HDF5_data_scheme") ; }   
   /// 
   bool has_key (std::string const & key) const { return (H5Lexists(this->getId(), key.c_str(), H5P_DEFAULT)); }
   ///
   void unlink_key_if_exists(std::string const & key) const { if (this->has_key(key)) this->unlink(key.c_str()); }
   /// Open a subgroup. Throw it if does not exists
   group open_group(std::string const & key) const { 
    if (!has_key(key))  TRIQS_RUNTIME_ERROR << "no subgroup "<<key <<" in the group";
    group res; 
    try { res = this->openGroup(key.c_str());} 
    catch (H5::GroupIException e){ TRIQS_RUNTIME_ERROR << "Error in opening the subgroup "<< key <<"\n H5 error message : \n "<< e.getCDetailMsg(); } 
    return res;
   }
   /// Open an existing DataSet. Throw it if does not exists
   H5::DataSet open_dataset(std::string const & key) const { 
    if (!has_key(key))  TRIQS_RUNTIME_ERROR << "no dataset "<<key <<" in the group";
    H5::DataSet res;
    try { res = this->openDataSet(key.c_str());} 
    catch (H5::GroupIException e){ TRIQS_RUNTIME_ERROR << "Error in opening the dataset "<< key <<"\n H5 error message : \n "<< e.getCDetailMsg(); } 
    return res;
   }
   /** 
    * \brief Create a subgroup.
    * \param key : the name of the subgroup
    * \param delete_if_exists : unlink the group if it exists 
    */
   group create_group(std::string const & key, bool delete_if_exists = true) const {
    unlink_key_if_exists(key); 
    return this->createGroup(key.c_str());
   }
   /** 
    * \brief Create a dataset.
    * \param key : the name of the subgroup
    * \param all others are forwarded to H5::Group
    * It unlinks the dataset if it exists.
    */
   template<typename ... Args>
    H5::DataSet create_dataset(std::string const & key, Args && ... args) const {
     unlink_key_if_exists(key); 
     H5::DataSet res;
     try { res = this->createDataSet(key.c_str(), std::forward<Args>(args)...);} 
     catch (H5::GroupIException e){ TRIQS_RUNTIME_ERROR << "Error in creating the dataset "<< key <<"\n H5 error message : \n "<< e.getCDetailMsg(); } 
     return res;
    }
 };
}}
#endif
