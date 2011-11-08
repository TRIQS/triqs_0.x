#  Copyright Olivier Parcollet 2010 
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

FILE(GLOB_RECURSE SOURCES *.rst)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in ${CMAKE_CURRENT_BINARY_DIR}/BUILD/conf.py)

SET(DOC_EXT "_cpp")
include(${TRIQS_SOURCE_DIR}/cmake/BuildSphinx.cmake)

add_subdirectory(mctools)
add_subdirectory(arrays)

# we must first have doxygen process the sources to link to them
add_dependencies(docs_sphinx_cpp docs_doxy_mctools)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/BUILD/html/ COMPONENT documentation DESTINATION share/doc/triqs/developer_manual
 FILES_MATCHING 
 PATTERN "*.html" 
 PATTERN "*.png"
 PATTERN "*.js"
 PATTERN "_*"
 PATTERN "*.png" 
 PATTERN "*.gif" 
 PATTERN "*.xsl" 
 PATTERN "*.css"
 PATTERN "*.pdf"
 PATTERN "*.py"
 PATTERN "*.txt"
 PATTERN "*.bib"
 PATTERN ".svn" EXCLUDE 
 PATTERN "CVS" EXCLUDE
 )

EXECUTE_PROCESS(COMMAND sh -c "mkdir -p BUILD && cd BUILD && ln -fs ${CMAKE_CURRENT_SOURCE_DIR}/* ${CMAKE_CURRENT_BINARY_DIR}/BUILD" )
