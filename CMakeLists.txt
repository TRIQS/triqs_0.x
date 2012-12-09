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
# and because I have not checked it (copy of various file before compilation may mess it up), blitz in particular....
IF (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
 MESSAGE(FATAL_ERROR "I am sorry, the on-site compilation is disabled at the moment. Use another directory")
ENDIF (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})

# nasty policy for escaping in string....
cmake_policy(SET CMP0011 OLD)
cmake_policy(SET CMP0005 NEW)
cmake_policy(SET CMP0003 NEW) #cmake --help-policy CMP0003 

cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
project (AllPrograms C CXX Fortran) 

# Make sure the rpath is set fine after installation
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Include TRIQS cmake macros
find_package(TriqsMacros)

# Check user input for BOOST_SOURCE_DIR
if(BOOST_SOURCE_DIR)
 if(NOT IS_ABSOLUTE ${BOOST_SOURCE_DIR})
  Message(FATAL_ERROR "You must specify an absolute path for BOOST_SOURCE_DIR.")
 endif(NOT IS_ABSOLUTE ${BOOST_SOURCE_DIR})
 STRING(REGEX MATCH "~\\/" TILDE "${BOOST_SOURCE_DIR}")
 IF(TILDE)
  Message(FATAL_ERROR "You should not use the tilde symbol ~. Please write the full absolute path for BOOST_SOURCE_DIR.")
 ENDIF(TILDE)
endif(BOOST_SOURCE_DIR)


# find the basic libraries from which TRIQS depends and add the def, include and link
find_package(TRIQS_dep_libs)
# Add macro definitions needed to use ALPS and dependent libraries
add_definitions(${ALPS_EXTRA_DEFINITIONS} ) 
# Add include directories needed to use ALPS and dependent libraries
include_directories(${ALPS_INCLUDE_DIRS} ${ALPS_EXTRA_INCLUDE_DIRS})
# Add link directories needed to use ALPS and dependent libraries
link_directories(${ALPS_LIBRARY_DIRS})
# reporting
Message(STATUS " Triqs general lib detection ( but boost ) ")
Message(STATUS " Adding definitions ${ALPS_EXTRA_DEFINITIONS}" )
Message(STATUS " Adding include  ${ALPS_INCLUDE_DIRS} ${ALPS_EXTRA_INCLUDE_DIRS} ")
Message(STATUS " Adding libs ${ALPS_EXTRA_LIBRARIES} " )
Message(STATUS " Adding library dir ${ALPS_LIBRARY_DIRS} " )

# Final destination of the python modules
string(REGEX REPLACE ".*/lib" "lib" TRIQS_PYTHON_LIB_DEST_ROOT ${ALPS_PYTHON_SITE_PKG} ) 
SET( TRIQS_PYTHON_LIB_DEST ${TRIQS_PYTHON_LIB_DEST_ROOT}/pytriqs ) 
MESSAGE(STATUS "Python modules will be installed in ${CMAKE_INSTALL_PREFIX}/${TRIQS_PYTHON_LIB_DEST}")

# compile boost if necessary
if(BOOST_SOURCE_DIR)
 add_subdirectory(foreignlibs/boost)
 MESSAGE(STATUS  "BOOST_INCLUDE_DIR = ${BOOST_INCLUDE_DIR}")
 include_directories(SYSTEM ${BOOST_INCLUDE_DIR})
else(BOOST_SOURCE_DIR)
 FIND_PACKAGE(Boost 1.46.0 REQUIRED )
 include_directories(${BOOST_INSTALL_ROOT_DIR}/include)
 set(BOOST_LIBRARY ${BOOST_PYTHON_LIB} ${BOOST_SERIALIZATION_LIB} ${BOOST_FILESYSTEM_LIB} ${BOOST_MPI_LIB})
endif(BOOST_SOURCE_DIR)

# find cblas.
find_package(CBLAS)
include_directories(SYSTEM ${CBLAS_INCLUDE_DIR})

Find_package(Git3) # our local version of FindGit...
if (GIT_FOUND)
 GIT_TREE_INFO ( ${CMAKE_SOURCE_DIR} "GIT")
 MESSAGE(STATUS "Git hash ${GIT_HASH}")
 MESSAGE(STATUS "Git branches : \n  ${GIT_BRANCHES}")
 add_definitions( -DGIT_HASH="${GIT_HASH}") # for all compilations
endif (GIT_FOUND)

# Compute WhoAmI and hostname
# Change this : make a version module
EXECUTE_PROCESS(COMMAND hostname OUTPUT_VARIABLE myHostName RESULT_VARIABLE returncode OUTPUT_STRIP_TRAILING_WHITESPACE) 
EXECUTE_PROCESS(COMMAND id -un OUTPUT_VARIABLE WhoAmI RESULT_VARIABLE returncode OUTPUT_STRIP_TRAILING_WHITESPACE)
add_definitions( 
 -D__COMPILEDBY__=${WhoAmI}
 -D__HOSTNAME__=${myHostName}
 )

# Adding 
add_definitions( -DBOOST_PP_VARIADICS=1 )

# Python support ?
#FIRST VERSION : if OFF, we should not try to detect the python at all....
# TEST AND DEBUG THE OFF OPTION
option(PythonSupport "Build with Python support (i.e. pytriqs, etc..) : OFF NOT TESTED ?" ON)
if (PythonSupport)
 add_definitions (-DTRIQS_WITH_PYTHON_SUPPORT)

 # Tool to help to include the wrappers in main and declare the built-in module
 Find_package(PythonModuleEmbeddingTool)
endif (PythonSupport)

enable_testing()

# on 64 bit machines
if (CMAKE_SIZEOF_VOID_P EQUAL 8) # for 64 bits machines
 MESSAGE( STATUS "64 bit machine : Adding -fpic ")
 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic" )
 add_definitions( -fpic)
endif (CMAKE_SIZEOF_VOID_P EQUAL 8)

#------------------------
# Documentation
##------------------------
option(Build_Documentation "Build the documentation ?" OFF)

# Built blitz and get corresponding variables
add_subdirectory(${TRIQS_SOURCE_DIR}/foreignlibs/blitz blitz)
add_definitions( -DUSE_BLITZ)

# remove the possible horrible pthread bug on os X !!
# check for clang compiler ?? on gcc, os X snow leopard, it MUST be set 
# since _REENTRANT is mysteriously set and blitz++ compiles in thread_safe mode
# and this leads to random stalling of the code....
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
 add_definitions( -pthread )
ENDIF( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# General include header
# remove this dep to C++
include_directories(${TRIQS_SOURCE_DIR})
include_directories(SYSTEM ${TRIQS_SOURCE_DIR}/foreignlibs ${BLITZ_INCLUDE} ${FFTW_INCLUDE_DIR})

# Create a temporary copy of the python modules so that we can run before installation with the test
#FILE(COPY ${TRIQS_SOURCE_DIR}/pytriqs DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ FILES_MATCHING PATTERN *.py PATTERN *.svn EXCLUDE )
if(BOOST_SOURCE_DIR)
 FILE(COPY ${TRIQS_SOURCE_DIR}/foreignlibs/boost DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ FILES_MATCHING PATTERN *.py PATTERN *.svn EXCLUDE )
else(BOOST_SOURCE_DIR)
 if(NOT BOOST_MODULE_DIR)
  EXEC_PYTHON_SCRIPT ("import boost; print boost.__path__[0]" BOOST_MODULE_DIR)
 endif(NOT BOOST_MODULE_DIR)
 execute_process(COMMAND ln -fs ${BOOST_MODULE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/pytriqs/boost)
 install (DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/pytriqs/boost DESTINATION ${TRIQS_PYTHON_LIB_DEST})
endif(BOOST_SOURCE_DIR)

#------------------------
# link libs for shared libs or for pytriqs 
#------------------------

SET(TRIQS_LINK_LIBS 
 ${BLITZ_LIBRARY}
 ${PYTHON_LIBRARY} ${PYTHON_EXTRA_LIBS} 
 ${FFTW_LIBRARIES}
 ${BOOST_LIBRARY} ${ALPS_EXTRA_LIBRARIES}
 ${LAPACK_LIBS}
 )

#------------------------
# General triqs libraries
#------------------------

add_subdirectory(${TRIQS_SOURCE_DIR}/triqs triqs)  

IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
 list (REMOVE_DUPLICATES TRIQS_LINK_LIBS)
ENDIF( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

#------------------------
# Applications
#------------------------

add_subdirectory(${TRIQS_SOURCE_DIR}/applications )

#------------------------
# pytriqs modules
#------------------------
if (PythonSupport)
 add_subdirectory(${TRIQS_SOURCE_DIR}/pytriqs ) 
endif (PythonSupport)

#------------------------
# build pytriqs main or make the script
#------------------------

add_subdirectory(shells bin)

#------------------------
# Documentation
#------------------------

#option(Build_Documentation "Build the documentation ?" OFF)
if (Build_Documentation)
 add_subdirectory (${TRIQS_SOURCE_DIR}/doc )
 #add_dependencies(docs_sphinx py_sources)
endif (Build_Documentation)

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

