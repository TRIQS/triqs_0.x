
SET(PYTHON_SOURCES 
  ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
  ${CMAKE_CURRENT_SOURCE_DIR}/gf_base.py
  ${CMAKE_CURRENT_SOURCE_DIR}/inverse.py
  ${CMAKE_CURRENT_SOURCE_DIR}/array_view.py
  ${CMAKE_CURRENT_SOURCE_DIR}/block_gf.py
  ${CMAKE_CURRENT_SOURCE_DIR}/gf_init.py
  ${CMAKE_CURRENT_SOURCE_DIR}/descriptors.py
  ${CMAKE_CURRENT_SOURCE_DIR}/lazy_expressions.py
  ${CMAKE_CURRENT_SOURCE_DIR}/gf_legendre.py
  ${CMAKE_CURRENT_SOURCE_DIR}/gf_imfreq.py
  ${CMAKE_CURRENT_SOURCE_DIR}/gf_imtime.py
  ${CMAKE_CURRENT_SOURCE_DIR}/gf_refreq.py
  ${CMAKE_CURRENT_SOURCE_DIR}/gf_retime.py
)

# Install python sources
install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/base/gf_local)

# Build python module
SET(SOURCES_WRAP gf_wrap.cpp)
python_build_module(GF _pytriqs_GF base/gf_local ${SOURCES_WRAP})

if (NOT TRIQS_BUILD_STATIC)
 target_link_libraries(_pytriqs_GF triqs ${TRIQS_LINK_LIBS})
endif (NOT TRIQS_BUILD_STATIC)


