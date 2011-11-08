#  Copyright Olivier Parcollet 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)


find_package(Doxygen REQUIRED)

# Prepare the Doxyfile
configure_file(${TRIQS_SOURCE_DIR}/cmake/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

SET(doxy_top ${CMAKE_CURRENT_BINARY_DIR}/xml_doxy)
#SET(DOXYGEN_EXECUTABLE "${DOXYGEN_EXECUTABLE} -u")
add_custom_command (OUTPUT ${doxy_top} DEPENDS ${DOXYGEN_SOURCES} COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
add_custom_target(docs_doxy ALL DEPENDS ${doxy_top})


