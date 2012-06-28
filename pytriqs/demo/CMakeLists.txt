
SET(PYTHON_SOURCES 
  ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
  ${CMAKE_CURRENT_SOURCE_DIR}/solve_python.py
)

install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/demo)


SET(CPP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/C++)

SET(SOURCES_CPP
 ${CMAKE_CURRENT_SOURCE_DIR}/mymodule_wrap.cpp
 )

python_build_module(demo_python mymodule demo ${SOURCES_CPP})



