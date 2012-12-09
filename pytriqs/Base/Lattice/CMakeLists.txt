
SET(PYTHON_SOURCES 
  ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
  ${CMAKE_CURRENT_SOURCE_DIR}/BZ_Patch.py
  ${CMAKE_CURRENT_SOURCE_DIR}/SuperLattice.py
  ${CMAKE_CURRENT_SOURCE_DIR}/TightBinding.py
)

install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/Base/Lattice)

# THE C++ code
SET(CPP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/C++)

SET(SOURCES_HPP 
 ${CPP_DIR}/bravais_lattice_and_brillouin_zone.hpp
 ${CPP_DIR}/functors.hpp
 ${CPP_DIR}/tight_binding.hpp
 )

SET(SOURCES_CPP
 ${CPP_DIR}/bravais_lattice_and_brillouin_zone.cpp
 ${CPP_DIR}/tight_binding.cpp
 ${CPP_DIR}/wrap.cpp
 )

python_build_module(LatticeTools _pytriqs_LatticeTools Base/Lattice ${SOURCES_CPP})



