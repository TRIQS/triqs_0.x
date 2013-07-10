# Path for local cmake modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Where is triqs source 
set(TRIQS_SOURCE_DIR "${CMAKE_SOURCE_DIR}")

# a default version : by default installation in the build dir : 
# useful for developer : everything is at the same place, easy to clean.
# pass -DCMAKE_INSTALL_PREFIX to overrule
set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/INSTALL_DIR" CACHE PATH "TRIQS install prefix" )
MESSAGE(STATUS "Installation directory is ${CMAKE_INSTALL_PREFIX}")

# Forbid compilation in the source dir... to force good practice
# and because I have not checked it (copy of various file before compilation may mess it up)
IF (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
 MESSAGE(FATAL_ERROR "I am sorry, the on-site compilation is disabled at the moment. Use another directory")
ENDIF (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})

# nasty policy for escaping in string....
cmake_policy(SET CMP0011 OLD)
cmake_policy(SET CMP0005 NEW)
cmake_policy(SET CMP0003 NEW) #cmake --help-policy CMP0003 

cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
project (AllPrograms C CXX Fortran) 

## C & C++ Headers
include(CheckIncludeFile)
include(CheckIncludeFileCXX)

# Make sure the rpath is set fine after installation
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# First check the C++ compiler
find_package(CompilerCheck)

# Include TRIQS cmake macros
find_package(TriqsMacros)

# Static built ? TO BE RECHECKED ...
option(BUILD_SHARED_LIBS "Build shared libraries" ON)
string(COMPARE EQUAL ${BUILD_SHARED_LIBS}  "OFF" TRIQS_BUILD_STATIC)
if (TRIQS_BUILD_STATIC)
 MESSAGE(STATUS "STATIC Built  ")
else (TRIQS_BUILD_STATIC) 
 MESSAGE(STATUS "DYNAMIC Built ")
endif (TRIQS_BUILD_STATIC) 

# RPATH setting (before compiling boost)
if(APPLE)
 set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
else(APPLE)
 set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
 set(CMAKE_SKIP_BUILD_RPATH FALSE)
 set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
 set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif(APPLE)

# #################################
# Python support 
# #################################
# TEST AND DEBUG THE OFF OPTION
option(PythonSupport "Build with Python support (i.e. pytriqs, etc..) : OFF NOT TESTED ?" ON)
if (PythonSupport)

 set(TRIQS_WITH_PYTHON_SUPPORT 1) # for the triqs_config.h file configuration

 find_package(Python REQUIRED)
 IF (NOT PYTHONLIBS_FOUND)
  MESSAGE (FATAL_ERROR "Can not find python")
 ENDIF()

 include_directories(SYSTEM ${PYTHON_NUMPY_INCLUDE_DIR} ${PYTHON_INCLUDE_DIRS})
 link_libraries( ${PYTHON_LIBRARY})

 MESSAGE (STATUS "Numpy include in ${PYTHON_NUMPY_INCLUDE_DIR}")

 # The C API of numpy has changed with 1.7.0, the macro is a version switch in a few files of the libs.
 if (PYTHON_NUMPY_VERSION VERSION_LESS "1.7.0")
  SET(TRIQS_NUMPY_VERSION_LT_17 1)
 endif()

 SET (TRIQS_LIBRARY_PYTHON ${PYTHON_LIBRARY} ${PYTHON_EXTRA_LIBS})
 # check : not a flag
 SET (TRIQS_INCLUDE_PYTHON ${PYTHON_INCLUDE_DIRS} ${PYTHON_NUMPY_INCLUDE_DIR})

 # Installation : Final destination of the python modules
 string(REGEX REPLACE ".*/lib" "lib" TRIQS_PYTHON_LIB_DEST_ROOT ${PYTHON_SITE_PKG} ) 
 SET( TRIQS_PYTHON_LIB_DEST ${TRIQS_PYTHON_LIB_DEST_ROOT}/pytriqs ) 
 MESSAGE(STATUS "Python modules will be installed in ${CMAKE_INSTALL_PREFIX}/${TRIQS_PYTHON_LIB_DEST}")

 # build pytriqs main or make the script
 add_subdirectory(shells bin)

endif (PythonSupport)

# #################################
# MPI
# #################################
MESSAGE( STATUS "-------- MPI detection -------------")
find_package(MPI)
if (NOT MPI_CXX_FOUND)
 MESSAGE(FATAL_ERROR "TRIQS requires MPI")
endif()
message (STATUS "MPI CXX compiler : ${MPI_CXX_COMPILER}")

MESSAGE( STATUS "MPI_COMPILE_FLAGS : ${MPI_CXX_COMPILE_FLAGS}")
MESSAGE( STATUS "MPI_CXX_INCLUDE_PATH : ${MPI_CXX_INCLUDE_PATH}")
add_definitions( ${MPI_CXX_COMPILE_FLAGS})
add_definitions( ${MPI_C_COMPILE_FLAGS})
include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})
link_libraries( ${MPI_CXX_LIBRARIES})
SET(TRIQS_LIBRARY_MPI ${MPI_CXX_LIBRARIES})
SET(TRIQS_INCLUDE_MPI ${MPI_CXX_INCLUDE_PATH})


# #################################
# Boost detection
# #################################

MESSAGE( STATUS "-------- Boost detection -------------")
# Check user input for BOOST_SOURCE_DIR
if(BOOST_SOURCE_DIR)
 if(NOT IS_ABSOLUTE ${BOOST_SOURCE_DIR})
  Message(FATAL_ERROR "You must specify an absolute path for BOOST_SOURCE_DIR.")
 endif(NOT IS_ABSOLUTE ${BOOST_SOURCE_DIR})
 STRING(REGEX MATCH "~\\/" TILDE "${BOOST_SOURCE_DIR}")
 IF(TILDE)
  Message(FATAL_ERROR "You should not use the tilde symbol ~. Please write the full absolute path for BOOST_SOURCE_DIR.")
 ENDIF(TILDE)

 # compile boost if necessary
 set(ALPS_BOOST_LIBRARY_NAME "boost_for_triqs" ) 
 SET(ALPS_HAVE_MPI 1)
 SET(ALPS_HAVE_BOOST_MPI 1)
 set(BUILD_BOOST_MPI TRUE)
 set(BUILD_BOOST_PYTHON TRUE)
 add_subdirectory(foreignlibs/boost)
 include_directories(SYSTEM ${BOOST_INCLUDE_DIR})

else(BOOST_SOURCE_DIR)

 if (NOT Boost_ROOT_DIR_found)
  SET( Boost_ROOT_DIR ${BOOST_ROOT_DIR})
 endif (NOT Boost_ROOT_DIR_found)

 # if static build, request static lib
 IF(TRIQS_BUILD_STATIC)
  set(Boost_USE_STATIC_LIBS        ON)
 else()
  set(Boost_USE_STATIC_LIBS        OFF)
 endif()
 # icc requires a higher version
 IF(CMAKE_COMPILER_IS_ICC)
  SET(BOOST_VERSION 1.53)
 else()
  SET(BOOST_VERSION 1.46)
 endif()
 if (PythonSupport)
  FIND_PACKAGE(Boost ${BOOST_VERSION} REQUIRED COMPONENTS mpi python serialization ) #filesystem  )
 Else()
  FIND_PACKAGE(Boost ${BOOST_VERSION} REQUIRED COMPONENTS mpi serialization ) #filesystem  )
 endif()
 if (NOT Boost_FOUND)
  MESSAGE(FATAL_ERROR "Boost not found")
 endif()
 MESSAGE(STATUS "BOOST has been detected")
 include_directories(${Boost_INCLUDE_DIRS})
 #include_directories(${BOOST_INSTALL_ROOT_DIR}/include)
 set(BOOST_LIBRARY ${Boost_LIBRARIES}) #${BOOST_PYTHON_LIB} ${BOOST_SERIALIZATION_LIB} ${BOOST_FILESYSTEM_LIB} ${BOOST_MPI_LIB})
 set(BOOST_INCLUDE_DIR ${Boost_INCLUDE_DIRS})
endif(BOOST_SOURCE_DIR)

MESSAGE(STATUS  "Boost include dir : ${BOOST_INCLUDE_DIR}")
MESSAGE(STATUS  "Boost libraries : ${BOOST_LIBRARY}")

#------------------------
# copy boost into pytriqs/boost
#------------------------
if (PythonSupport)
 if(NOT BOOST_SOURCE_DIR)
  if(NOT BOOST_MODULE_DIR)
   EXEC_PYTHON_SCRIPT ("import boost; print boost.__path__[0]" BOOST_MODULE_DIR)
  endif(NOT BOOST_MODULE_DIR)
  execute_process(COMMAND ln -fs ${BOOST_MODULE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/pytriqs/)
  install (DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/pytriqs/boost DESTINATION ${TRIQS_PYTHON_LIB_DEST})
 endif(NOT BOOST_SOURCE_DIR)

endif (PythonSupport)


# on 64 bit machines
if (CMAKE_SIZEOF_VOID_P EQUAL 8) # for 64 bits machines
 MESSAGE( STATUS "64 bit machine : Adding -fpic ")
 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic" )
 add_definitions( -fpic)
endif (CMAKE_SIZEOF_VOID_P EQUAL 8)

# #################################
# Lapack
# #################################
MESSAGE( STATUS "-------- Lapack detection -------------")
#option(LAPACK_64_BIT "Use 64-bit version of LAPACK" OFF)
#set (BIND_FORTRAN_INTEGER_8 ${LAPACK_64_BIT})
if (NOT LAPACK_FOUND)
 find_package(Lapack)
endif (NOT LAPACK_FOUND)
IF (REQUIRE_PTHREAD)
 SET(LAPACK_LIBRARY "${LAPACK_LIBRARY};${PTHREAD_LIBRARY}")
 SET(LAPACK_LIBRARIES "${LAPACK_LIBRARIES};${PTHREAD_LIBRARY}")
ENDIF (REQUIRE_PTHREAD)

SET(LAPACK_LIBS ${LAPACK_LIBRARY} ${BLAS_LIBRARY} ${LAPACK_LINKER_FLAGS} CACHE STRING "Flags to link Lapack and Blas")

SET(TRIQS_LIBRARY_LAPACK ${LAPACK_LIBS})
if(MAC_VECLIB)
    set(TRIQS_LIBRARY_LAPACK "-framework vecLib ${TRIQS_LIBRARY_LAPACK}")
endif(MAC_VECLIB) 
MESSAGE(STATUS "Lapack libraries : ${TRIQS_LIBRARY_LAPACK}")
link_libraries(${TRIQS_LIBRARY_LAPACK})

# #################################
# HDF5
# #################################
MESSAGE( STATUS "-------- HDF5 detection -------------")
# on weiss, it is 2.8.2 and we should not put HL, on 12.04 we need to put it...
if ( ${CMAKE_VERSION} VERSION_LESS "2.8.6") # CHECK THIS BOUND, where are the cmake changelogs ??
 find_package(HDF5 REQUIRED C CXX )
else(${CMAKE_VERSION} VERSION_LESS "2.8.6")
 find_package(HDF5 REQUIRED C CXX HL )
endif(${CMAKE_VERSION} VERSION_LESS "2.8.6")
IF(NOT HDF5_FOUND)
 MESSAGE(FATAL_ERROR "Require hdf5 1.8.2 or higher. Set HDF5_HOME")
ENDIF()
IF(HDF5_IS_PARALLEL)
 MESSAGE(FATAL "parallel(MPI) hdf5 is detected. The standard version is preferred.")
ENDIF(HDF5_IS_PARALLEL)
add_definitions (${HDF5_DEFINITIONS})
include_directories (${HDF5_INCLUDE_DIR})
link_libraries (${HDF5_LIBRARIES} ) #${HDF5_CXX_LIBRARIES} )
MESSAGE( STATUS " HDF5_LIBRARIES = ${HDF5_LIBRARIES} ")

# #################################
# Other libs
# #################################
# GMP
find_package(GMP)

# FFTW
find_package(FFTW)


#
## Remove the test from the build ?
option(Build_Triqs_General_Tools_Test "Build the test of the general tool ? ( array class , ... ) " ON)

# TRIQS cmake file and config.h
get_property(ALL_DEFS GLOBAL PROPERTY TRIQS_COMPIL_DEFINITIONS)
configure_file(./cmake/triqs_config.h.in Config/triqs_config.h @ONLY)
IF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
 configure_file(./cmake/TRIQSConfig.cmake.in Config/TRIQSConfig.cmake @ONLY)
else()
 configure_file(./cmake/TRIQSConfig.cmake.in Config/TRIQSConfig1.cmake @ONLY )
 EXEC_PYTHON_SCRIPT( "import re;open('${CMAKE_CURRENT_BINARY_DIR}/Config/TRIQSConfig.cmake','w').write(re.sub(r';*\\s*-framework\\s+([a-zA-Z_][a-zA-Z_0-9]*)',r' \"-framework \\1\"',open('${CMAKE_CURRENT_BINARY_DIR}/Config/TRIQSConfig1.cmake').read()))", ConfigFile)
endif()

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/Config/TRIQSConfig.cmake DESTINATION share/triqs/cmake)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/Config/triqs_config.h DESTINATION include/triqs)

# 64 bits machine ?
#if(CMAKE_SIZEOF_VOID_P EQUAL 8 OR CMAKE_GENERATOR MATCHES Win64)
#  set (ALPS_64BIT ON)
#endif(CMAKE_SIZEOF_VOID_P EQUAL 8 OR CMAKE_GENERATOR MATCHES Win64)

# is it useful ?
if(NOT CMAKE_BUILD_TYPE)
 set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
 mark_as_advanced(CMAKE_BUILD_TYPE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

# Fortran bindings
option(BIND_FORTRAN_LOWERCASE "FORTRAN functions are compiled WITHOUT a trailing underscore" OFF)
mark_as_advanced(BIND_FORTRAN_LOWERCASE)
if (BIND_FORTRAN_LOWERCASE)
 SET(TRIQS_BIND_FORTRAN_LOWERCASE OFF )
endif (BIND_FORTRAN_LOWERCASE)

# git
Find_package(Git3) # our local version of FindGit...
if (GIT_FOUND)
 GIT_TREE_INFO ( ${CMAKE_SOURCE_DIR} "GIT")
 MESSAGE(STATUS "Git hash ${GIT_HASH}")
 MESSAGE(STATUS "Git branches : \n  ${GIT_BRANCHES}")
endif (GIT_FOUND)

# Compute WhoAmI and hostname
# Change this : make a version module
EXECUTE_PROCESS(COMMAND hostname OUTPUT_VARIABLE myHostName RESULT_VARIABLE returncode OUTPUT_STRIP_TRAILING_WHITESPACE) 
EXECUTE_PROCESS(COMMAND id -un OUTPUT_VARIABLE WhoAmI RESULT_VARIABLE returncode OUTPUT_STRIP_TRAILING_WHITESPACE)
SET(TRIQS_HOSTNAME ${myHostName})
SET(TRIQS_COMPILEDBY ${WhoAmI})

#
enable_testing()

##------------------------
# Documentation
##------------------------
option(Build_Documentation "Build the documentation ?" OFF)
option(Build_PDF_Documentation "Build the pdf documentation ?" OFF)

# remove the possible horrible pthread bug on os X !!
# check for clang compiler ?? on gcc, os X snow leopard, it MUST be set 
# since _REENTRANT is mysteriously set and this leads to random stalling of the code....
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
 add_definitions( -pthread )
ENDIF( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# General include header
# remove this dep to C++
include_directories(${TRIQS_SOURCE_DIR})
include_directories(SYSTEM ${FFTW_INCLUDE_DIR})

#------------------------
# link libs for shared libs or for pytriqs 
#------------------------

SET(TRIQS_LINK_LIBS 
 ${TRIQS_LIBRARY_PYTHON}
 ${FFTW_LIBRARIES}
 ${BOOST_LIBRARY} 
 ${LAPACK_LIBS}
 ${GMP_LIBRARIES} ${GMPXX_LIBRARIES}
 )

#------------------------
# General triqs libraries
#------------------------

add_definitions ( -DTRIQS_BUILDING_LIBRARY -I ${CMAKE_BINARY_DIR}/Config)
add_subdirectory(${TRIQS_SOURCE_DIR}/triqs triqs)  

IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
 list (REMOVE_DUPLICATES TRIQS_LINK_LIBS)
ENDIF( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

#------------------------
# Applications
#------------------------

#add_subdirectory(${TRIQS_SOURCE_DIR}/applications )

#------------------------
# pytriqs modules
#------------------------
if (PythonSupport)
 add_subdirectory(${TRIQS_SOURCE_DIR}/pytriqs ) 
endif (PythonSupport)

#------------------------
# tests
#------------------------
add_subdirectory(${TRIQS_SOURCE_DIR}/test ) 

#------------------------
# Documentation
#------------------------

#option(Build_Documentation "Build the documentation ?" OFF)
if (Build_Documentation)
 add_subdirectory (${TRIQS_SOURCE_DIR}/doc )
 #add_dependencies(docs_sphinx py_sources)
endif (Build_Documentation)

#------------------------
# FindTRIQS
#------------------------

install (FILES ${CMAKE_SOURCE_DIR}/cmake/FindTRIQS.cmake DESTINATION share/triqs/cmake)

#------------------------
# Packaging
#------------------------

# make the source package
SET(TAR_SRC_FILE ${CMAKE_CURRENT_BINARY_DIR}/triqs_src.tgz)
SET(TAR_SRC_DIRBUILD ${CMAKE_CURRENT_BINARY_DIR}/triqs_src_tgz_build)
add_custom_command( OUTPUT ${TAR_SRC_FILE} COMMAND mkdir -p  ${TAR_SRC_DIRBUILD} && cd ${TAR_SRC_DIRBUILD} && ${Subversion_SVN_EXECUTABLE} export --force ${Project_WC_URL} triqs_src && tar cvzf ${TAR_SRC_FILE} triqs_src)
add_custom_target(src_tgz DEPENDS ${TAR_SRC_FILE})

# .deb and dmg
INCLUDE(InstallRequiredSystemLibraries)

# Prepare the postinstall and preremoval for debian
configure_file(${TRIQS_SOURCE_DIR}/cmake/prerm.in ${CMAKE_CURRENT_BINARY_DIR}/prerm)
configure_file(${TRIQS_SOURCE_DIR}/cmake/postinst.in ${CMAKE_CURRENT_BINARY_DIR}/postinst)

SET(CPACK_PACKAGE_NAME "triqs")
SET(CPACK_PACKAGE_FILE_NAME "triqs")
SET(CPACK_PACKAGE_VENDOR "The TRIQS collaboration")
SET(CPACK_DEBIAN_PACKAGE_SECTION "Science")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "TRIQS, a Toolbox for Research in Quantum Interacting Systems")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.txt")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")
SET(CPACK_PACKAGE_CONTACT "Olivier Parcollet <olivier.parcollet@cea.fr>")
SET(CPACK_PACKAGE_VERSION_MAJOR "0")
SET(CPACK_PACKAGE_VERSION_MINOR "9")
SET(CPACK_PACKAGE_VERSION_PATCH "1-2")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "TRIQS_${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_STRIP_FILES "bin/pytriqs")
SET(CPACK_SOURCE_STRIP_FILES "")
SET(CPACK_PACKAGE_EXECUTABLES "pytriqs" "pytriqs")
SET(CPACK_DEBIAN_PACKAGE_DEPENDS "openmpi-bin(>=1.4), libblas-dev (>=1.2), ipython (>=0.12), libfftw3-dev (>=3.2), hdf5-tools (>=1.8), python-matplotlib (>=1.1) ")
SET(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA ${CMAKE_CURRENT_BINARY_DIR}/prerm;${CMAKE_CURRENT_BINARY_DIR}/postinst)
IF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
 SET(CPACK_GENERATOR TGZ DEB)
ENDIF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

INCLUDE(CPack)

