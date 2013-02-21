#include "./parameters.hpp"

namespace triqs { namespace utility {
 std::map<size_t, std::function<_object(std::string const &)>> _object::deserialization_fnts;
 std::map<size_t, std::function<_object(group_or_file const &, std::string const &)>> _object::h5_read_fnts;
 bool parameters::initialized = false;
}}

