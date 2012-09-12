# list the sources
SET(SOURCES_HPP 
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/Hloc.hpp  
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/SmallMatrix.hpp  
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/TraceSlice.hpp  
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/TraceSliceStack.hpp  
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/TimeEvolution.hpp  
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/Time_Ordered_Operator_List.hpp  
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/DynamicTrace.hpp  
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/Configuration.hpp  
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/detManip.hpp  
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/det_manip.hpp  
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/gf_binner_and_eval.hpp  
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/Move_AuxiliaryFunctions.hpp  
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/Move_Global.hpp  
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/Move_Insert_Remove_Cdag_C_Delta.hpp  
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/Move_Insert_Remove_Cdag_C_Delta_SegmentPicture.hpp  
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/Move_Move_CC_Delta.hpp  
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/myConstIteratorVector.hpp  
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/util.hpp
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/Matsubara_generators.hpp  
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/Measures_F.hpp  
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/Measures_G.hpp  
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/Measures_Legendre.hpp  
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/Measures_Legendre_allseries.hpp  
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/Measures_OpAv.hpp  
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/Measures_OpCorr.hpp  
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/Measures_Z.hpp  
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/MC.hpp  
 )

SET(SOURCES_CPP 
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/Hloc.cpp 
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/DynamicTrace.cpp 
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/Measures_OpCorr.cpp   
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/Measures_Legendre.cpp 
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/Measures_Legendre_allseries.cpp 
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/Configuration.cpp 
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/MC.cpp 
 ${CMAKE_CURRENT_SOURCE_DIR}/C++/QMC_wrap.cpp
 )

SET(PYTHON_SOURCES 
 ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/Solver.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/Models.py 
 )

#SET(CODENAME cthyb)
#include(${CMAKE_SOURCE_DIR}/cmake/PrintEditCode.cmake)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/C++)

option( Build_CTHyb "Should I build CTHyb ?" ON) 
if (Build_CTHyb)
 python_build_module(CTHyb "_pytriqs_Solver_HybridizationExpansion" Solvers/HybridizationExpansion ${SOURCES_CPP})
 install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/Solvers/HybridizationExpansion)
 add_subdirectory(test)
endif (Build_CTHyb)

#set(NickName "CTHyb")
#python_build_module_optional (${NickName}  "_pytriqs_Solver_HybridizationExpansion" Solvers/HybridizationExpansion ${SOURCES_CPP} )
#if (Build_${NickName})
# if (NOT TRIQS_BUILD_STATIC) 
#  target_link_libraries (_pytriqs_Solver_HybridizationExpansion triqs)
# endif (NOT TRIQS_BUILD_STATIC) 
# install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/Solvers/HybridizationExpansion)
# add_subdirectory(test)
#endif (Build_${NickName})

