#  Copyright Olivier Parcollet 2010 
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

SET(SOURCES ${SOURCES}
${CMAKE_CURRENT_SOURCE_DIR}/advanced.rst
${CMAKE_CURRENT_SOURCE_DIR}/interface.rst
${CMAKE_CURRENT_SOURCE_DIR}/LDADMFTmain.rst
${CMAKE_CURRENT_SOURCE_DIR}/selfcons.rst
${CMAKE_CURRENT_SOURCE_DIR}/analysis.rst
${CMAKE_CURRENT_SOURCE_DIR}/wien2k.rst
${CMAKE_CURRENT_SOURCE_DIR}/Ce-HI.rst
PARENT_SCOPE
)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/BUILD/html COMPONENT documentation DESTINATION share/doc/triqs/wien2k
 FILES_MATCHING 
 PATTERN "*.pdf" 
)

