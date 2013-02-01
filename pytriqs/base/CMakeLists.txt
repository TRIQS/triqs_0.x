
configure_file ( ${CMAKE_CURRENT_SOURCE_DIR}/version.py.in ${CMAKE_CURRENT_BINARY_DIR}/version.py)

SET(PYTHON_SOURCES 
  ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
  ${CMAKE_CURRENT_BINARY_DIR}/version.py
  )

install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/base)

add_subdirectory(archive)
add_subdirectory(dmft)
add_subdirectory(dos)
add_subdirectory(gf)
add_subdirectory(gf_local)
add_subdirectory(plot)
add_subdirectory(fit)
add_subdirectory(utility)

add_subdirectory(lattice)
add_subdirectory(sumk)

add_subdirectory(test)

