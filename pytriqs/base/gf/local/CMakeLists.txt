
SET(PYTHON_SOURCES 
  ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
)

# Install python sources
install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/base/gf)

# Build python module
#SET(SOURCES_WRAP gf_wrap.cpp)
#python_build_module(GF2 _pytriqs_GF2 base/gf/local ${SOURCES_WRAP})

#if (NOT TRIQS_BUILD_STATIC)
# target_link_libraries(_pytriqs_GF2 triqs ${TRIQS_LINK_LIBS})
#endif (NOT TRIQS_BUILD_STATIC)

# Build cython module 

cython_module( GF3 gf base/gf/local )



