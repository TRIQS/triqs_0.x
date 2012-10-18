
SET(PYTHON_SOURCES 
  ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
  ${CMAKE_CURRENT_SOURCE_DIR}/TailGF.py
)

# Install python sources
install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/Base/gf)

# Build python module
#SET(SOURCES_WRAP gf_wrap.cpp)
#python_build_module(GF2 _pytriqs_GF2 Base/gf/local ${SOURCES_WRAP})

#if (NOT TRIQS_BUILD_STATIC)
# target_link_libraries(_pytriqs_GF2 triqs ${TRIQS_LINK_LIBS})
#endif (NOT TRIQS_BUILD_STATIC)

# Build cython module 

#cython_module( GF3 _pytriqs_GF3 Base/gf/local ${CMAKE_CURRENT_SOURCE_DIR}/TailGF.pyx )



