SET(ExecutableName "pytriqs" ) 

if (TRIQS_BUILD_STATIC)

 MESSAGE( STATUS "Preparing the executables")
 # The main target. 
 add_executable(${ExecutableName} ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
 install(TARGETS ${ExecutableName} DESTINATION bin)  
 include_directories( ${CMAKE_CURRENT_BINARY_DIR} )

 # Find if ipython is present
 find_package(IPython)
 if (IPYTHON_AVAIL) 
  MESSAGE(STATUS "ipython is available on your platform: I will compile pytriqs and ipytriqs")
  add_executable( i${ExecutableName} ${CMAKE_CURRENT_SOURCE_DIR}/main_i.cpp)
  install(TARGETS i${ExecutableName} DESTINATION bin)  
 endif (IPYTHON_AVAIL) 

 # incorporate the correct pythonpath 
 SET(PythonPathAdd "${CMAKE_INSTALL_PREFIX}/${TRIQS_PYTHON_LIB_DEST_ROOT}")
 add_definitions( -DPYTHONPATH_ADD=${PythonPathAdd})

 # exploit the result of Modules about the lib to link and which module was compiled (because of options).
 GET_PROPERTY(DEPENDANCE_TO_ADD  GLOBAL  PROPERTY DEPENDANCE_TO_ADD)
 foreach (d ${DEPENDANCE_TO_ADD})
  add_dependencies(${ExecutableName} ${d} )
  if (IPYTHON_AVAIL) 
   add_dependencies(i${ExecutableName} ${d} )
  endif (IPYTHON_AVAIL) 
 endforeach (d ${DEPENDANCE_TO_ADD})

 # what should I link ? 
 GET_PROPERTY(LIBS_TO_LINK  GLOBAL  PROPERTY LIBS_TO_LINK)
 target_link_libraries( ${ExecutableName} ${LIBS_TO_LINK} ${TRIQS_LINK_LIBS} triqs)
 if (IPYTHON_AVAIL) 
  target_link_libraries( i${ExecutableName} ${LIBS_TO_LINK} ${TRIQS_LINK_LIBS} triqs)
 endif (IPYTHON_AVAIL) 

 # which modules did I build in ?
 GET_PROPERTY(Modules_actually_compiled  GLOBAL  PROPERTY Modules_actually_compiled)
 STRING( REPLACE ";" " " ModuleList "${Modules_actually_compiled}") 
 add_definitions( -DMODULE_LIST=${ModuleList})

 # generate the header to build the modules in 
 if(BOOST_SOURCE_DIR)
   set_property (GLOBAL APPEND PROPERTY PYTHON_STATIC_MODULES_LIST mpi)
 endif(BOOST_SOURCE_DIR)
 generate_python_header_for_main()

else (TRIQS_BUILD_STATIC)

 MESSAGE( STATUS "Preparing the various scripts ")
 # a simple script to set the python path and have the same effect as the statically compiled code.
 option(TRIQS_ADD_PYPATHS_TO_SHELLS "Add the correct python paths to shell (ON except for packaging)" ON)
 mark_as_advanced( TRIQS_ADD_PYPATHS_TO_SHELLS)
 if (TRIQS_ADD_PYPATHS_TO_SHELLS)
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/pytriqs_shell.bash.in ${CMAKE_CURRENT_BINARY_DIR}/pytriqs)
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/ipytriqs_shell.bash.in ${CMAKE_CURRENT_BINARY_DIR}/ipytriqs)
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/ipytriqs_shell.notebook.bash.in ${CMAKE_CURRENT_BINARY_DIR}/ipytriqs_notebook)
 else (TRIQS_ADD_PYPATHS_TO_SHELLS)
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/pytriqs_shell.nopypath.bash.in ${CMAKE_CURRENT_BINARY_DIR}/pytriqs)
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/ipytriqs_shell.nopypath.bash.in ${CMAKE_CURRENT_BINARY_DIR}/ipytriqs)
 endif (TRIQS_ADD_PYPATHS_TO_SHELLS)

 configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/pytriqs_for_test_shell.bash.in ${CMAKE_BINARY_DIR}/run_pytriqs_for_test.sh)
 execute_process( COMMAND chmod u+x  ${CMAKE_BINARY_DIR}/run_pytriqs_for_test.sh )

 install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pytriqs ${CMAKE_CURRENT_BINARY_DIR}/ipytriqs DESTINATION bin  
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE )

endif (TRIQS_BUILD_STATIC)


