
SET(PYTHON_SOURCES 
  ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
  ${CMAKE_CURRENT_SOURCE_DIR}/_GFBloc_concept_impl.py
  ${CMAKE_CURRENT_SOURCE_DIR}/_GFBloc_base_data_tail.py
  ${CMAKE_CURRENT_SOURCE_DIR}/inverse.py
  ${CMAKE_CURRENT_SOURCE_DIR}/ArrayViewWithIndexConverter.py
  ${CMAKE_CURRENT_SOURCE_DIR}/GF.py
  ${CMAKE_CURRENT_SOURCE_DIR}/GF_Initializers.py
  ${CMAKE_CURRENT_SOURCE_DIR}/Descriptors.py
  ${CMAKE_CURRENT_SOURCE_DIR}/lazy_expressions.py
  ${CMAKE_CURRENT_SOURCE_DIR}/GFBloc_ImLegendre.py
  ${CMAKE_CURRENT_SOURCE_DIR}/GFBloc_ImFreq.py
  ${CMAKE_CURRENT_SOURCE_DIR}/GFBloc_ImTime.py
  ${CMAKE_CURRENT_SOURCE_DIR}/GFBloc_ReFreq.py
  ${CMAKE_CURRENT_SOURCE_DIR}/GFBloc_ReTime.py
)

# Install python sources
install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/base/GF_Local)

# Build python module
SET(SOURCES_WRAP GF_wrap.cpp)
python_build_module(GF _pytriqs_GF base/GF_Local ${SOURCES_WRAP})

if (NOT TRIQS_BUILD_STATIC)
 target_link_libraries(_pytriqs_GF triqs ${TRIQS_LINK_LIBS})
endif (NOT TRIQS_BUILD_STATIC)


