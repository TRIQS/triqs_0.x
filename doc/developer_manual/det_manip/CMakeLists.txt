#  Copyright Olivier Parcollet 2010 
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

SET(SOURCES ${SOURCES}
${CMAKE_CURRENT_SOURCE_DIR}/det_manip.rst
PARENT_SCOPE
)

SET(DOXYGEN_SOURCES ${TRIQS_SOURCE_DIR}/triqs/det_manip/det_manip.hpp)
SET(DOXYGEN_NAME det_manip)
SET(DOXYGEN_HTML_OUTPUT ../html/doxy_detmanip)

include(${TRIQS_SOURCE_DIR}/cmake/BuildDoxygenDocs.cmake)

