#include "./parameters.hpp"

namespace triqs { namespace utility {

 std::map<size_t, std::function<_object(std::string const &)>> _object::deserialization_fnts;
 std::map<size_t, std::function<_object(group_or_file const &, std::string const &)>> _object::h5_read_fnts;
 std::map<hid_t, std::pair<size_t,std::string> > _object::h5_type_to_c_equivalence;
 bool _object::initialized = false;

 //-----------------------------------------------------------------------

 void _object::_init() { 
   if (initialized) return;
#define REGISTER_UNSERIALIZE_FUN(r, data, T) _object::register_native_type<T>(); 
   BOOST_PP_SEQ_FOR_EACH(REGISTER_UNSERIALIZE_FUN, nil , TRIQS_UTIL_PARAM_PREDEFINED_CAST); 
#undef REGISTER_UNSERIALIZE_FUN
   register_type<parameters>(); 
   initialized = true;
  }

 //-----------------------------------------------------------------------

  void h5_read ( group_or_file f, std::string const & name, _object & obj){
  using namespace H5; 
  if (!f.exists(name)) TRIQS_RUNTIME_ERROR << "no such "<<name <<" in file ";
  size_t type_hash;
  // First if it is not a subgroup, it is a scalar or a an array
  // TO BE CHECKED
  try {
   H5::DataSet ds = f->openDataSet( name.c_str() );
   hid_t dtype = H5Dget_type(ds.getId()); 
   hid_t native_type=H5Tget_native_type(dtype, H5T_DIR_DEFAULT);

/*   std::cerr  << H5Tequal(native_type, H5T_NATIVE_DOUBLE)<<std::endl ;
   std::cerr  << H5Tequal(native_type, H5T_NATIVE_LONG)<<std::endl ;
   std::cerr  << H5Tequal(native_type, H5T_NATIVE_INT)<<std::endl ;
   std::cerr  << H5Tequal(native_type,         H5T_NATIVE_CHAR         )<<std::endl ;
   std::cerr  << H5Tequal(dtype, H5T_NATIVE_INT)<<std::endl ;
   std::cerr  << H5Tequal(dtype, H5T_NATIVE_LONG)<<std::endl ;
   std::cerr  << H5Tequal(dtype, H5T_NATIVE_INT)<<std::endl ;
   std::cerr  << H5Tequal(dtype, H5T_C_S1)<<std::endl ;
   std::cerr<< (H5Tget_class(dtype) == H5T_STRING )<< std::endl ; 
*/

   if (H5Tget_class(dtype) == H5T_STRING ) { // it is a string
    type_hash = type_hash_code<std::string>(); 
   }
   else if ((H5Tget_class(dtype) == H5T_INTEGER ) || ((H5Tget_class(dtype) == H5T_FLOAT ) )) { 
    auto it = _object::h5_type_to_c_equivalence.begin();
    for (;it != _object::h5_type_to_c_equivalence.end();++it) { if (H5Tequal(native_type, it->first))  break;}
    if (it == _object::h5_type_to_c_equivalence.end()) TRIQS_RUNTIME_ERROR << " h5_type_to_c_equivalence : type not found";
    type_hash = it->second.first;
   }

   herr_t status = H5Tclose (native_type);
   status = H5Tclose (dtype);
  }
  TRIQS_ARRAYS_H5_CATCH_EXCEPTION;

  //to try only. Need to make h5_code function for all object....
  // write, read an attribute with the type ? 
  obj = _object::h5_read_fnts[type_hash](f,name);
  }

  //-----------------------------------------------------------------------

  // move this is the subgroup name
  extern "C" { 
   herr_t get_group_elements_name(hid_t loc_id, const char *name, void *opdata) { static_cast<std::vector<std::string> *>(opdata)->push_back(name); return 0;}
  }

  void h5_read ( group_or_file F, std::string const & subgroup_name, parameters & p){ 
   auto gr = F.open_group(subgroup_name);
   _object obj;
   std::vector<std::string> s; 
   //Use iterator to see the names of the objects in the file
   int ret = F->iterateElems(subgroup_name, NULL, get_group_elements_name, static_cast<void*>(&s));

   for (auto & x : s) {
    std::cerr  << " loading : "<< x <<std::endl;
    try {
     h5_read(gr,x,obj);
     p[x] = obj;
    }
    catch(...) { std::cerr  << " Can not load "<< x<<" !!"<< std::endl ;}
   }
  }

}}

