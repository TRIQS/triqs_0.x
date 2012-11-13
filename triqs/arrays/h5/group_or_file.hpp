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
#ifndef TRIQS_ARRAYS_H5_GROUP_OR_FILE_H
#define TRIQS_ARRAYS_H5_GROUP_OR_FILE_H
#include "./common.hpp"
namespace triqs { namespace arrays { namespace h5 { 

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

  group_or_file(std::string filename, int flag) { 
   H5::H5File f(filename, H5F_ACC_TRUNC);
   fg_ptr = boost::make_shared<H5::H5File>(f); id = f.getId();
  }

  group_or_file(hid_t id_) {
   H5::Group g; g.setId(id_); 
   fg_ptr = boost::make_shared<H5::Group>(g);
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

}}}
#endif
