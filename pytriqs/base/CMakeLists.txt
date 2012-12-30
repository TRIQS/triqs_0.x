
configure_file ( ${CMAKE_CURRENT_SOURCE_DIR}/version.py.in ${CMAKE_CURRENT_BINARY_DIR}/version.py)

SET(PYTHON_SOURCES 
  ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
  ${CMAKE_CURRENT_BINARY_DIR}/version.py
  )

install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/base)

add_subdirectory(Archive)
add_subdirectory(DMFT)
add_subdirectory(DOS)

option(BUILD_NEW_GF " Build the new GF (require the github version of cython ...)" OFF)
if (BUILD_NEW_GF) 
 add_subdirectory(gf)
endif(BUILD_NEW_GF)

add_subdirectory(GF_Local)
add_subdirectory(Plot)
add_subdirectory(Fit)
add_subdirectory(Utility)

add_subdirectory(Lattice)
add_subdirectory(SumK)

add_subdirectory(test)

