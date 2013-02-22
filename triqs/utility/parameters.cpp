#include "./parameters.hpp"

namespace triqs { namespace utility {

 std::map<size_t, std::function<_object(std::string const &)>> _object::deserialization_fnts;
 std::map<size_t, std::function<_object(h5::group const &, std::string const &)>> _object::h5_read_fnts;
 std::map<hid_t, std::pair<size_t,std::string> > _object::h5_type_to_c_equivalence;
 std::map<std::pair<size_t,int>, size_t > _object::hash_code_element_rank_2_hash_code_array;
 std::map<std::string, size_t > _object::h5_scheme_to_code;
 bool _object::initialized = false;

 //-----------------------------------------------------------------------

 void _object::_init() {
  if (initialized) return;
#define REGISTER_UNSERIALIZE_FUN(r, data, T) _object::register_native_type<T>();
  BOOST_PP_SEQ_FOR_EACH(REGISTER_UNSERIALIZE_FUN, nil , TRIQS_UTIL_PARAM_PREDEFINED_CAST);
#undef REGISTER_UNSERIALIZE_FUN
  initialized = true;
 }

 //-----------------------------------------------------------------------

 void h5_read ( h5::group g, std::string const & name, _object & obj){
  //std::cerr  << " hread object "<< name << std::endl ;
  using namespace H5;
  if (!g.has_key(name)) TRIQS_RUNTIME_ERROR << "no such "<<name <<" in file ";
  size_t type_hash;
  // First if it is not a subgroup, it is a scalar or a an array
  // TO BE CHECKED
  try {
   H5::DataSet ds = g.open_dataset( name.c_str() );
   hid_t dtype = H5Dget_type(ds.getId());
   hid_t native_type=H5Tget_native_type(dtype, H5T_DIR_DEFAULT);
   if (H5Tget_class(dtype) == H5T_STRING ) { // it is a string
    type_hash = type_hash_code<std::string>();
   }
   else if ((H5Tget_class(dtype) == H5T_INTEGER ) || ((H5Tget_class(dtype) == H5T_FLOAT ) )) {
    int rank = ds.getSpace().getSimpleExtentNdims();
    auto it = _object::h5_type_to_c_equivalence.begin();
    for (;it != _object::h5_type_to_c_equivalence.end();++it) { if (H5Tequal(native_type, it->first))  break;}
    if (it == _object::h5_type_to_c_equivalence.end()) TRIQS_RUNTIME_ERROR << " h5_type_to_c_equivalence : type not found";
    type_hash = it->second.first;
    if (rank>0) {
     size_t type_hash_element = type_hash;
     auto it2= _object::hash_code_element_rank_2_hash_code_array.find(std::make_pair(type_hash_element,rank));
     if (it2 == _object::hash_code_element_rank_2_hash_code_array.end())
      TRIQS_RUNTIME_ERROR << " hash_code_element_rank_2_hash_code_array : type not found" << rank;
     type_hash = it2->second;
    }
   }
   //else if (H5Tget_class(dtype)==H5T_ARRAY)
   herr_t status = H5Tclose (native_type);
   status = H5Tclose (dtype);
  }
  TRIQS_ARRAYS_H5_CATCH_EXCEPTION;

  obj = _object::h5_read_fnts[type_hash](g,name);
 }

 //-----------------------------------------------------------------------

 // move this is the subgroup name
 extern "C" {
  herr_t get_group_elements_name_ds(hid_t loc_id, const char *name, void *opdata) {
   H5O_info_t object_info; herr_t err =H5Oget_info_by_name(loc_id,name,&object_info,H5P_DEFAULT );
   if (err<0) TRIQS_RUNTIME_ERROR << "get_group_elements_name_ds internal";
   if (object_info.type == H5O_TYPE_DATASET)  static_cast<std::vector<std::string> *>(opdata)->push_back(name);
   return 0;
  }
  herr_t get_group_elements_name_grp(hid_t loc_id, const char *name, void *opdata) {
   H5O_info_t object_info; herr_t err =H5Oget_info_by_name(loc_id,name,&object_info,H5P_DEFAULT );
   if (err<0) TRIQS_RUNTIME_ERROR << "get_group_elements_name_grp internal";
   if (object_info.type == H5O_TYPE_GROUP)  static_cast<std::vector<std::string> *>(opdata)->push_back(name);
   return 0;
  }
 }
 //-----------------------------------------------------------------------

 void h5_read ( h5::group F, std::string const & subgroup_name, parameters & p){
  auto gr = F.open_group(subgroup_name);
  std::vector<std::string> ds_name, grp_name;
  //Use iterator to see the names of the objects in the file
  int ret = F.iterateElems(subgroup_name.c_str(), NULL, get_group_elements_name_ds, static_cast<void*>(&ds_name));
  ret = F.iterateElems(subgroup_name.c_str(), NULL, get_group_elements_name_grp, static_cast<void*>(&grp_name));
  for (auto & x : grp_name) {
   //std::cerr  << " loading group : "<< x <<std::endl;
   auto x_grp = gr.open_group(x);
   auto triqs_data_scheme = x_grp.read_triqs_hdf5_data_scheme(); 
   if (triqs_data_scheme != "") { 
    auto type_hash = _object::h5_scheme_to_code[triqs_data_scheme];
    auto it = _object::h5_read_fnts.find(type_hash);
    if (it == _object::h5_read_fnts.end()) TRIQS_RUNTIME_ERROR << "TRIQS_HDF5_data_scheme : ["<< triqs_data_scheme << "] is unknown. Did you register your object ?";
    p[x] = it->second(gr,x);
   }
   else { 
    parameters p2;
    h5_read (gr,x,p2);
    p[x] = p2;
   }
  }
  for (auto & x : ds_name) {
   //std::cerr  << " loading : "<< x <<std::endl;
   //try {
   _object obj;
   h5_read(gr,x,obj);
   p[x] = obj;
   // }
   //catch(...) { std::cerr  << " Can not load "<< x<<" !!"<< std::endl ;}
  }
 }

 //-----------------------------------------------------------------------

 void parameters::update (parameters const & pdef){ for (auto const & pvp : pdef) (*this)[pvp.first] = pvp.second; }

 void parameters::update_with_defaults (parameters const & pdef, ull_t flag ){

  if ( (flag & reject_key_without_default) ) { // check that no other parameter is present
   for (auto const & pvp : *this) 
    if (!pdef.has_key( pvp.first)) 
     TRIQS_RUNTIME_ERROR << "update : parameter "<< pvp.first << " is absent from the defaults and no_parameter_without_default is ON. ";
  }

  for (auto const & pvp : pdef) {
   if (this->has_key( pvp.first)) { // check the type is correct 
    if (pvp.second.type_num != (*this)[pvp.first].type_num)
     TRIQS_RUNTIME_ERROR << "update_with_defaults : parameter "<< pvp.first << " does not have the right type";
   }
   else { 
    (*this)[pvp.first] = pvp.second; // insert the default
   }
  }
 }

}}

