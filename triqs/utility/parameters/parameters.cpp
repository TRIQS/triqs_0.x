#include "./parameters.hpp"

namespace triqs { namespace utility {

 void h5_read ( h5::group F, std::string const & subgroup_name, parameters & p){
  auto gr = F.open_group(subgroup_name);
  std::vector<std::string> ds_name = F.get_all_dataset_names(subgroup_name), grp_name = F.get_all_subgroup_names(subgroup_name);
  for (auto & x : grp_name) {
   //std::cerr  << " loading group : "<< x <<std::endl;
   auto x_grp = gr.open_group(x);
   auto triqs_data_scheme = x_grp.read_triqs_hdf5_data_scheme(); 
   if (triqs_data_scheme != "") { 
    auto type_hash = _object::h5_scheme_to_code[triqs_data_scheme];
    auto it = _object::code_to_h5_read_fnts.find(type_hash);
    if (it == _object::code_to_h5_read_fnts.end()) TRIQS_RUNTIME_ERROR << "TRIQS_HDF5_data_scheme : ["<< triqs_data_scheme << "] is unknown. Did you register your object ?";
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
   try {
    _object obj;
    h5_read(gr,x,obj);
    p[x] = obj;
   }
   catch(H5::Exception const & e) { TRIQS_RUNTIME_ERROR<< "Can not load "<< x<<"\n H5 error is : \n   "<< e.getCDetailMsg();}
  }
 }

 //-----------------------------------------------------------------------

 void parameters::update (parameters const & pdef){ for (auto const & pvp : pdef) (*this)[pvp.first] = pvp.second; }

 void parameters::update (parameters_defaults const & pdef, ull_t flag ){

  if ( (flag & reject_key_without_default) ) { // check that no other parameter is present
   for (auto const & pvp : *this) 
    if (!pdef.has_key( pvp.first)) 
     TRIQS_RUNTIME_ERROR << "update : parameter "<< pvp.first << " is absent from the defaults and no_parameter_without_default is ON. ";
  }

  for (auto const & pvp : pdef) {
   if (this->has_key( pvp.first)) { // check the type is correct 
    if (pvp.second.type_code_ != (*this)[pvp.first].type_code_)
     TRIQS_RUNTIME_ERROR << "update_with_defaults : parameter "<< pvp.first << " does not have the right type";
   }
   else { 
    (*this)[pvp.first] = pvp.second; // insert the default
   }
  }
 }

}}

