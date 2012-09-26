#  Copyright Olivier Parcollet 2010 
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

# Sphinx sources
#FILE(GLOB_RECURSE SOURCES *.rst)
#set_property(GLOBAL APPEND PROPERTY SPHINX_SOURCES ${SOURCES}) 

# Doxygen sources
set_property(GLOBAL APPEND PROPERTY DOXYGEN_SOURCES 
 ${TRIQS_SOURCE_DIR}/triqs/arrays/h5/group_or_file.hpp
 ${TRIQS_SOURCE_DIR}/triqs/arrays/h5/simple_read_write.hpp
 ${TRIQS_SOURCE_DIR}/triqs/arrays/h5/array_proxy.hpp
 ${TRIQS_SOURCE_DIR}/triqs/arrays/h5/array_stack.hpp
 ${TRIQS_SOURCE_DIR}/triqs/arrays/array.hpp
 )

