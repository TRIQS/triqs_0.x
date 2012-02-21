
find_package(Doxygen REQUIRED)

string(REPLACE ";" "  " DOXYGEN_SOURCES_LIST "${DOXYGEN_SOURCES}")
#message( ${DOXYGEN_SOURCES}) 

# Prepare the Doxyfile
configure_file(${TRIQS_SOURCE_DIR}/cmake/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

SET(doxy_top ${CMAKE_CURRENT_BINARY_DIR}/doxy.log)
add_custom_command (OUTPUT ${doxy_top} DEPENDS ${DOXYGEN_SOURCES} COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile > ${doxy_top} )
add_custom_target(docs_doxy_${DOXYGEN_NAME} ALL DEPENDS ${doxy_top})

