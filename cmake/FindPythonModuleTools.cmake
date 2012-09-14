#  Copyright Olivier Parcollet 2012 
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

# Tools to make building python module easy

# --------- C++ ------------------

# Compile and link with python 
#link_libraries(${PYTHON_LIBRARY} ${PYTHON_EXTRA_LIBS} )

# Why So General ?? It is included in all codes !!!
include_directories(${PYTHON_INCLUDE_DIRS} ${PYTHON_NUMPY_INCLUDE_DIR})

# This function add the target to build a python module
#
# NickName = 
# ModuleName = 
# ModuleDest = path in the pytriqs tree [ FOR INSTALLATION ONLY] 
#
function (python_build_module NickName ModuleName ModuleDest  )
 MESSAGE(STATUS "Preparing python module  ${NickName} ")
 
 # prepare a simple .py file, the same name for static and dynamic case....
 STRING(REGEX REPLACE "^_" "" MODNAME ${ModuleName})
 SET(ModuleIncludeFile pytriqs/${ModuleDest}/${MODNAME}.py)
 MESSAGE (STATUS "Making ${ModuleIncludeFile}")

 if (TRIQS_BUILD_STATIC)
 
  add_library(triqs_${NickName} ${ARGN} )
  target_link_libraries (${ModuleName} triqs)

  # VIRER Modules_actually_compiled : useless now... 
  set_property (GLOBAL APPEND PROPERTY Modules_actually_compiled ${NickName} )
  set_property (GLOBAL APPEND PROPERTY PYTHON_STATIC_MODULES_LIST ${ModuleName})
  set_property (GLOBAL APPEND PROPERTY LIBS_TO_LINK triqs_${NickName} )
  file (WRITE ${CMAKE_BINARY_DIR}/${ModuleIncludeFile}  "from _${MODNAME} import *")
 
 else (TRIQS_BUILD_STATIC)

  set_property (GLOBAL APPEND PROPERTY PYTHON_DYNAMIC_MODULES_LIST ${ModuleName})
  add_library(${ModuleName} MODULE ${ARGN}  )
  set_target_properties(${ModuleName}  PROPERTIES PREFIX "") #eliminate the lib in front of the module name 
  target_link_libraries(${ModuleName} ${TRIQS_LINK_LIBS} triqs)
  install (TARGETS ${ModuleName} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/${ModuleDest}  )

  # SHould be in the static part ??? Useless ??
  set_property (GLOBAL APPEND PROPERTY DEPENDANCE_TO_ADD triqs_${NickName} )

  STRING(REPLACE "/" "." MODPATH ${ModuleDest})
  # issue #26: setting RTLD_GLOBAL is no longer necessary because boost python is loaded dynamically and it was causing problems
  # with a double definition of PyArray_API from multiarray.so
  # file (WRITE ${CMAKE_BINARY_DIR}/${ModuleIncludeFile}  "${FIX}\nimport sys,ctypes \nflag = sys.getdlopenflags()\nsys.setdlopenflags(flag|ctypes.RTLD_GLOBAL)\nfrom pytriqs.${MODPATH}._${MODNAME} import *\nsys.setdlopenflags(flag)")
  file (WRITE ${CMAKE_BINARY_DIR}/${ModuleIncludeFile}  "from pytriqs.${MODPATH}._${MODNAME} import *\n")
 
 endif (TRIQS_BUILD_STATIC)
 
 install( FILES ${CMAKE_BINARY_DIR}/${ModuleIncludeFile} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/${ModuleDest})
endfunction (python_build_module ModuleName)


# REPLACE AND DELETE THIS ....

#function (python_build_module_optional_OBSOL NickName ModuleName ModuleDest)
# MESSAGE(STATUS "Preparing optional python module  ${NickName} ")
# option( Build_${NickName} "Should I build ${NickName} ?" ON) #OFF)
# if (Build_${NickName})
#  python_build_module(${NickName} ${ModuleName} ${ModuleDest} ${ARGN}) 
# endif (Build_${NickName})
#endfunction (python_build_module_optional ModuleNickName ModuleName ModuleDest)

# Only used in f90 ?
# Move the OPTION OUT of HERE !
function (python_register_dynamic_module NickName)
 option( Build_${NickName} "Should I build ${NickName} ?" ON) #OFF)
 if (Build_${NickName})
  set_property (GLOBAL APPEND PROPERTY DEPENDANCE_TO_ADD triqs_${NickName} )
  set_property (GLOBAL APPEND PROPERTY Modules_actually_compiled ${NickName} )
 endif (Build_${NickName})
endfunction (python_register_dynamic_module NickName)

# Only for static case...
if (TRIQS_BUILD_STATIC)
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
endif (TRIQS_BUILD_STATIC)

#----------- F90 -------------------

#
# This macro builds the f2py module
#   - target_name
#   - 
#
macro (BuildF2pyModule target_name modulename module_pyf_name filelist1)

 set(filelist ${filelist1}  ${ARGN})
 set(filename temp_script.py)
 # Copy all the files
 EXECUTE_PROCESS(COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/${module_pyf_name} ${CMAKE_CURRENT_BINARY_DIR} )
 FOREACH( f ${filelist})
  EXECUTE_PROCESS(COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/${f} ${CMAKE_CURRENT_BINARY_DIR} )
 ENDFOREACH(f)
 # write the script that will build the f2py extension
 SET(filename ${CMAKE_CURRENT_BINARY_DIR}/${filename} )
 FILE(WRITE ${filename} "import sys\n")
 FILE(APPEND ${filename} "from numpy.f2py import main\n")
 FILE(APPEND ${filename} "sys.argv = [''] +'-c --f90exec=${CMAKE_Fortran_COMPILER} -m ${modulename} ${modulename}.pyf ${filelist} -llapack'.split()\n") 
 FILE(APPEND ${filename} "main()\n")

 # We had the normal target of the module
 add_custom_target(${target_name} ALL DEPENDS ${modulename}.so) 

 # TODO : to be corrected with the filelist is more than one file.
 # ... and a special target to build vertex.so, that depends on the sources files
 add_custom_command(OUTPUT  ${modulename}.so 
  COMMAND echo See `pwd`/f2pyBuild.log for logs
  COMMAND ${PYTHON_INTERPRETER} temp_script.py > f2pyBuild.log 2>&1
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${filelist} ${CMAKE_CURRENT_SOURCE_DIR}/${module_pyf_name} 
  )

 #FILE(RELATIVE_PATH rel ${CMAKE_SOURCE_DIR}/Modules ${CMAKE_CURRENT_SOURCE_DIR}/)
 #install (FILES ${CMAKE_CURRENT_BINARY_DIR}/${modulename}.so
 # DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/${ExecutableName}/triqs/${rel}/..)

endmacro (BuildF2pyModule)



