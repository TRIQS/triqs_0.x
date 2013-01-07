set(NickName "DFT")

# install the python sources 
SET(PYTHON_SOURCES 
 ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/solver_multiband.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/sumk_lda.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/sumk_lda_tools.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/symmetry.py 
 ${CMAKE_CURRENT_SOURCE_DIR}/U_matrix.py 
 )

SET( D ${CMAKE_CURRENT_SOURCE_DIR}/dmftproj/SRC_templates/)
SET(WIEN_SRC_TEMPL_FILES ${D}/case.cf_f_mm2  ${D}/case.cf_p_cubic  ${D}/case.indmftpr ${D}/run_triqs ${D}/runsp_triqs)

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
 install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/dft)
 install (FILES ${WIEN_SRC_TEMPL_FILES} DESTINATION share/triqs/Wien2k_SRC_files/SRC_templates )
 add_subdirectory(dmftproj) 
 add_subdirectory(F90)
 add_subdirectory(test)
 add_subdirectory(converters)
endif (Build_${NickName})
