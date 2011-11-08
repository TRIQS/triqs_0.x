
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

#
#find_package(InstalltriqsModule)
#InstalltriqsModule()
install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/Base/GF_Local)

# THE C++ code
SET(CPP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/C++)

SET(SOURCES_HPP 
 ${CPP_DIR}/TailGF.hpp 
 ${CPP_DIR}/MeshGF.hpp 
 ${CPP_DIR}/fourier.hpp 
 ${CPP_DIR}/legendre_matsubara.hpp 
 ${CPP_DIR}/GF_Bloc_Base.hpp
 ${CPP_DIR}/GF_Bloc_ImFreq.hpp
 ${CPP_DIR}/GF_Bloc_ReFreq.hpp
 ${CPP_DIR}/GF_Bloc_ImTime.hpp 
 ${CPP_DIR}/GF_Bloc_ReTime.hpp 
 ${CPP_DIR}/GF_Bloc_ImLegendre.hpp 
 )

SET(SOURCES_CPP 
 ${CPP_DIR}/TailGF.cpp 
 ${CPP_DIR}/MeshGF.cpp 
 ${CPP_DIR}/fourier.cpp 
 ${CPP_DIR}/legendre_matsubara.cpp 
 ${CPP_DIR}/GF_Bloc_Base.cpp
 ${CPP_DIR}/GF_Bloc_ImFreq.cpp
 ${CPP_DIR}/GF_Bloc_ReFreq.cpp
 ${CPP_DIR}/GF_Bloc_ImTime.cpp 
 ${CPP_DIR}/GF_Bloc_ReTime.cpp 
 ${CPP_DIR}/GF_Bloc_ImLegendre.cpp 
)

SET(SOURCES_WRAP
 ${CPP_DIR}/GF_wrap.cpp
 )

SET(CODENAME gf)
include(${CMAKE_SOURCE_DIR}/cmake/PrintEditCode.cmake)

add_library(LocalGreenFunction_c  ${SOURCES_CPP}  )
python_build_module(GF _pytriqs_GF Base/GF_Local ${SOURCES_WRAP})

if (NOT TRIQS_BUILD_STATIC)
 #add_library(LocalGreenFunction_c SHARED ${SOURCES_CPP}  )
 target_link_libraries(LocalGreenFunction_c ${TRIQS_LINK_LIBS})
 target_link_libraries(_pytriqs_GF LocalGreenFunction_c ${TRIQS_LINK_LIBS})
endif (NOT TRIQS_BUILD_STATIC)


