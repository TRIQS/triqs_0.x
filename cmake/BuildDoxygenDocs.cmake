
find_package(Doxygen REQUIRED)

string(REPLACE ";" "  " DOXYGEN_SOURCES_LIST "${DOXYGEN_SOURCES}")
#message( ${DOXYGEN_SOURCES}) 

# Prepare the Doxyfile
SET(DOXY_MCTOOLS_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/../BUILD/html/doxy_${DOXYGEN_NAME})
configure_file(${TRIQS_SOURCE_DIR}/cmake/Doxyfile.in ${DOXY_MCTOOLS_BUILD_DIR}/Doxyfile)

SET(doxy_top ${DOXY_MCTOOLS_BUILD_DIR}/doxy.log)
add_custom_command (OUTPUT ${doxy_top} DEPENDS ${DOXYGEN_SOURCES} COMMAND cd ${DOXY_MCTOOLS_BUILD_DIR} && ${DOXYGEN_EXECUTABLE} Doxyfile > ${doxy_top} )
add_custom_target(docs_doxy_${DOXYGEN_NAME} ALL DEPENDS ${doxy_top})



