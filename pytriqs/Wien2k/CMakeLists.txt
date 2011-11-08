set(NickName "Wien2k")

# install the python sources 
SET(PYTHON_SOURCES 
 ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/Solver_MultiBand.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/SumK_LDA_SO.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/SumK_LDA_SOtools.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/SumK_LDA.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/SumK_LDA_tools.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/Symmetry.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/Symmetry_SO.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/SumK_LDA_Wien2k_input.py
 ${CMAKE_CURRENT_SOURCE_DIR}/Umatrix.py 
 )

SET( D ${CMAKE_CURRENT_SOURCE_DIR}/dmftproj/SRC_templates/)
SET(WIEN_SRC_TEMPL_FILES ${D}/case.cf_f_mm2  ${D}/case.cf_p_cubic  ${D}/case.indmftpr)

#
python_register_dynamic_module(${NickName})

# build the fortran stuff...
if (Build_${NickName})
 message(STATUS "-----------------------------------------------------------------------------")
 message(STATUS "                             ********  WARNING ******** ")
 message(STATUS "  Wien2k users : after installation of TRIQS, copy the files from ")
 message(STATUS "  ${CMAKE_INSTALL_PREFIX}/share/triqs/Wien2k_SRC_files/SRC_templates ")
 message(STATUS "  to your Wien2k installation WIENROOT/SRC_templates (Cf documentation).     ")
 message(STATUS "  This is not handled automatically by the installation process. ")
 message(STATUS "-----------------------------------------------------------------------------")
 install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/Wien2k)
 install (FILES ${WIEN_SRC_TEMPL_FILES} DESTINATION share/triqs/Wien2k_SRC_files/SRC_templates )
 add_subdirectory(dmftproj) 
 add_subdirectory(F90)
 add_subdirectory(test)
endif (Build_${NickName})
