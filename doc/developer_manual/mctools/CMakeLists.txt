#  Copyright Olivier Parcollet 2010 
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

# Sphinx sources
#set_property(GLOBAL APPEND PROPERTY SPHINX_SOURCES  
#${CMAKE_CURRENT_SOURCE_DIR}/contents.rst
#${CMAKE_CURRENT_SOURCE_DIR}/introduction.rst
#${CMAKE_CURRENT_SOURCE_DIR}/move_concept.rst
#${CMAKE_CURRENT_SOURCE_DIR}/measure_concept.rst
#${CMAKE_CURRENT_SOURCE_DIR}/example.rst
#${CMAKE_CURRENT_SOURCE_DIR}/mc_generic.rst
#)

# Doxygen sources
set_property(GLOBAL APPEND PROPERTY DOXYGEN_SOURCES ${TRIQS_SOURCE_DIR}/triqs/mc_tools/mc_generic.hpp)
