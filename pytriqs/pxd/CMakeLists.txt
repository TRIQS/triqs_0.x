FILE(GLOB all_pxd_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.pxd )
install (FILES ${all_pxd_files} DESTINATION "include/pytriqs")
