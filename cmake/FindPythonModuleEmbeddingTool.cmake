
# Compile and link with python 
#link_libraries(${PYTHON_LIBRARY} ${PYTHON_EXTRA_LIBS} )
include_directories(${PYTHON_INCLUDE_DIRS} ${PYTHON_NUMPY_INCLUDE_DIR})

function (python_build_module NickName ModuleName ModuleDest  )
MESSAGE(STATUS "Preparing python module  ${NickName} ")
set_property (GLOBAL APPEND PROPERTY Modules_actually_compiled ${NickName} )
STRING(REGEX REPLACE "^_" "" MODNAME ${ModuleName})
SET(ModuleIncludeFile pytriqs/${ModuleDest}/${MODNAME}.py)
MESSAGE (STATUS "Making ${ModuleIncludeFile}")
if (TRIQS_BUILD_STATIC)
 add_library(triqs_${NickName} ${ARGN} )
 set_property (GLOBAL APPEND PROPERTY PYTHON_STATIC_MODULES_LIST ${ModuleName})
 set_property (GLOBAL APPEND PROPERTY LIBS_TO_LINK triqs_${NickName} )
 file (WRITE ${CMAKE_BINARY_DIR}/${ModuleIncludeFile}  "from _${MODNAME} import *")
else (TRIQS_BUILD_STATIC)
 set_property (GLOBAL APPEND PROPERTY PYTHON_DYNAMIC_MODULES_LIST ${ModuleName})
 add_library(${ModuleName} MODULE ${ARGN}  )
 set_target_properties(${ModuleName}  PROPERTIES PREFIX "") #eliminate the lib in front of the module name 
 target_link_libraries(${ModuleName} ${TRIQS_LINK_LIBS} triqs)
 install (TARGETS ${ModuleName} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/${ModuleDest}  )
 set_property (GLOBAL APPEND PROPERTY DEPENDANCE_TO_ADD triqs_${NickName} )
 STRING(REPLACE "/" "." MODPATH ${ModuleDest})
  # issue #26: setting RTLD_GLOBAL is no longer necessary because boost python is loaded dynamically and it was causing problems
  # with a double definition of PyArray_API from multiarray.so
  # file (WRITE ${CMAKE_BINARY_DIR}/${ModuleIncludeFile}  "${FIX}\nimport sys,ctypes \nflag = sys.getdlopenflags()\nsys.setdlopenflags(flag|ctypes.RTLD_GLOBAL)\nfrom pytriqs.${MODPATH}._${MODNAME} import *\nsys.setdlopenflags(flag)")
  file (WRITE ${CMAKE_BINARY_DIR}/${ModuleIncludeFile}  "from pytriqs.${MODPATH}._${MODNAME} import *\n")
endif (TRIQS_BUILD_STATIC)
install( FILES ${CMAKE_BINARY_DIR}/${ModuleIncludeFile} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/${ModuleDest})
endfunction (python_build_module ModuleName)

function (python_build_module_optional NickName ModuleName ModuleDest)
MESSAGE(STATUS "Preparing optional python module  ${NickName} ")
option( Build_${NickName} "Should I build ${NickName} ?" ON) #OFF)
if (Build_${NickName})
 python_build_module(${NickName} ${ModuleName} ${ModuleDest} ${ARGN}) 
endif (Build_${NickName})
endfunction (python_build_module_optional ModuleNickName ModuleName ModuleDest)

function (python_register_dynamic_module NickName)
option( Build_${NickName} "Should I build ${NickName} ?" ON) #OFF)
if (Build_${NickName})
 set_property (GLOBAL APPEND PROPERTY DEPENDANCE_TO_ADD triqs_${NickName} )
 set_property (GLOBAL APPEND PROPERTY Modules_actually_compiled ${NickName} )
endif (Build_${NickName})
endfunction (python_register_dynamic_module NickName)

function (generate_python_header_for_main)
GET_PROPERTY(PYTHON_STATIC_MODULES_LIST  GLOBAL  PROPERTY PYTHON_STATIC_MODULES_LIST)
# The files to be included in the main for the wrapping.
SET (Main_includewrap_file    ${CMAKE_CURRENT_BINARY_DIR}/main_includewrap.cpp)
SET (Main_appendinittab_file  ${CMAKE_CURRENT_BINARY_DIR}/main_appendinittab.cpp)
FILE( REMOVE ${Main_includewrap_file} ${Main_appendinittab_file}) 
foreach (module_name ${PYTHON_STATIC_MODULES_LIST})
 file (APPEND ${Main_includewrap_file}   "extern \"C\" { void init${module_name}();}; \n")
 file (APPEND ${Main_appendinittab_file} "PyImport_AppendInittab((char*)(\"${module_name}\"), &init${module_name});\n")
endforeach (module_name ${PYTHON_STATIC_MODULES_LIST})
endfunction (generate_python_header_for_main)


# This function add the target to build a python module
#
# NickName =  
# ModuleName = the python name of the module 
# ModuleDest = path in the pytriqs tree [ FOR INSTALLATION ONLY] 
macro (cython_module NickName ModuleName ModuleDest  )
 MESSAGE(STATUS "Preparing cython module  ${NickName} ")
 get_filename_component(CYTHON_EXECUTABLE_PATH ${PYTHON_INTERPRETER} PATH)
 SET(CYTHON_EXECUTABLE ${CYTHON_EXECUTABLE_PATH}/cython CACHE STRING "Path to the cython executable")
 SET(cython_src ${CMAKE_CURRENT_SOURCE_DIR}/${ModuleName}.pyx )
 FILE(GLOB all_pyx_src RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.pyx *.pxd )
 SET(cython_wrap  ${CMAKE_CURRENT_BINARY_DIR}/wrap_${NickName}_by_cython.cpp)
 add_custom_command (OUTPUT ${cython_wrap} DEPENDS ${all_pyx_src} ${ARGN} COMMAND ${CYTHON_EXECUTABLE} ${cython_src} -I ${CMAKE_CURRENT_SOURCE_DIR}/ --cplus -o ${cython_wrap}  )
 add_custom_target(cython_${NickName} ALL DEPENDS ${cython_wrap})
 python_build_module(${NickName} ${ModuleName} ${ModuleDest} ${cython_wrap} )
endmacro (cython_module)

