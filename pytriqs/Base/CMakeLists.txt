
configure_file ( ${CMAKE_CURRENT_SOURCE_DIR}/version.py.in ${CMAKE_CURRENT_BINARY_DIR}/version.py)

SET(PYTHON_SOURCES 
  ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
  ${CMAKE_CURRENT_BINARY_DIR}/version.py
  )

install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/Base)

add_subdirectory(Archive)
add_subdirectory(DMFT)
add_subdirectory(DOS)
add_subdirectory(GF_Local)
add_subdirectory(Plot)
add_subdirectory(Fit)
add_subdirectory(Utility)

add_subdirectory(Lattice)
add_subdirectory(SumK)

add_subdirectory(test)

#deprecated
#add_subdirectory (Lattice1)
#add_subdirectory(SumK1)

