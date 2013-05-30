#  Copyright Olivier Parcollet 2011
#  Adapted from the alps cmakelist
#
#  Copyright Synge Todo and Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

# strange lower case compare to other options... changing this
if (NOT Boost_ROOT_DIR_found)
 SET( Boost_ROOT_DIR ${BOOST_ROOT_DIR})
endif (NOT Boost_ROOT_DIR_found)

# choose static or dynamic
option(BUILD_SHARED_LIBS "Build shared libraries" ON)
string(COMPARE EQUAL ${BUILD_SHARED_LIBS}  "OFF" TRIQS_BUILD_STATIC)
if (TRIQS_BUILD_STATIC)
 MESSAGE(STATUS "STATIC Built  ")
else (TRIQS_BUILD_STATIC) 
 MESSAGE(STATUS "DYNAMIC Built ")
endif (TRIQS_BUILD_STATIC) 

option(LAPACK_64_BIT "Use 64-bit version of LAPACK" OFF)
option(BIND_FORTRAN_LOWERCASE "FORTRAN functions are compiled without a trailing underscore" OFF)
mark_as_advanced(BIND_FORTRAN_LOWERCASE)

message (STATUS "CMAKE_GENERATOR: ${CMAKE_GENERATOR}")
message (STATUS "CMAKE_CL_64: ${CMAKE_CL_64}")
message (STATUS "CMAKE_SIZEOF_VOID_P: ${CMAKE_SIZEOF_VOID_P}")

if(CMAKE_SIZEOF_VOID_P EQUAL 8 OR CMAKE_GENERATOR MATCHES Win64)
  set (ALPS_64BIT ON)
endif(CMAKE_SIZEOF_VOID_P EQUAL 8 OR CMAKE_GENERATOR MATCHES Win64)

message (STATUS "CMAKE_GENERATOR: ${CMAKE_GENERATOR}")
message (STATUS "CMAKE_CL_64: ${CMAKE_CL_64}")

set(ALPS_BOOST_LIBRARY_NAME "boost_for_triqs" ) 
set(ALPS_INSTALL_HEADERS  ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
  mark_as_advanced(CMAKE_BUILD_TYPE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})


######################################################################
# C & C++ Headers
######################################################################
include(CheckIncludeFile)
include(CheckIncludeFileCXX)

######################################################################
# Packages
######################################################################

# Boost Libraries if existing. don't need that at the moment
#find_package(BoostForALPS REQUIRED)

# MPI
set(ALPS_ENABLE_MPI ON)
if(ALPS_ENABLE_MPI)
  IF(NOT MPI_FOUND)
    IF(ALPS_PACKAGE_LIBRARIES AND NOT ALPS_USE_VISTRAILS AND APPLE)
      set(MPI_COMPILER ${CMAKE_INSTALL_PREFIX}/bin/mpicxx)
      message (STATUS "MPI compiler is ${MPI_COMPILER}")
    ENDIF(ALPS_PACKAGE_LIBRARIES AND NOT ALPS_USE_VISTRAILS AND APPLE)
    find_package(MPI)
      message (STATUS "MPI compiler was ${MPI_COMPILER}")
  ENDIF(NOT MPI_FOUND)

  set(MPI_DEFINITIONS ${MPI_COMPILE_FLAGS})
  if(MPI_INCLUDE_DIR)
    set(MPI_DEFINITIONS "${MPI_COMPILE_FLAGS} -I${MPI_INCLUDE_DIR}")
  endif(MPI_INCLUDE_DIR)
 
  set(MPI_INCLUDE_DIR ${MPI_INCLUDE_PATH})
  IF(MPI_FOUND)
   include_directories(SYSTEM ${MPI_INCLUDE_DIR})
    SET(ALPS_HAVE_MPI 1)
    SET(ALPS_HAVE_BOOST_MPI 1)
    set(BUILD_BOOST_MPI TRUE)
    foreach(arg ${MPI_LIBRARIES})
      set(MPI_LIBS "${MPI_LIBS} ${arg}")
    endforeach(arg ${MPI_LIBRARIES})
  ELSE(MPI_FOUND)
    set(BUILD_BOOST_MPI FALSE)
  ENDIF(MPI_FOUND)
else (ALPS_ENABLE_MPI)
    SET(ALPS_HAVE_MPI 0)
    SET(ALPS_HAVE_BOOST_MPI 0)
    set(BUILD_BOOST_MPI FALSE)
endif(ALPS_ENABLE_MPI)

# OpenMP
option(ALPS_ENABLE_OPENMP "Enable OpenMP parallelization" OFF)
option(ALPS_ENABLE_OPENMP_WORKER "Enable OpenMP worker support" OFF)
mark_as_advanced(ALPS_ENABLE_OPENMP)
mark_as_advanced(ALPS_ENABLE_OPENMP_WORKER)
if(ALPS_ENABLE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    # Almost always OpenMP flags are same both for C and for Fortran.
    if(ALPS_BUILD_FORTRAN)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
    endif(ALPS_BUILD_FORTRAN)
    if(ALPS_ENABLE_OPENMP_WORKER)
      # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DALPS_ENABLE_OPENMP_WORKER")
    endif(ALPS_ENABLE_OPENMP_WORKER)
  endif(OPENMP_FOUND)
endif(ALPS_ENABLE_OPENMP)

# Lapack
if (NOT LAPACK_FOUND)
 find_package(Lapack)
endif (NOT LAPACK_FOUND)
if (LAPACK_FOUND)
 set(ALPS_HAVE_LAPACK 1)
 set(ALPS_HAVE_BLAS 1)
endif(LAPACK_FOUND)
if(MSVC)
 file(COPY ${LAPACK_LIBRARY} ${BLAS_LIBRARY} DESTINATION ${PROJECT_BINARY_DIR}/lib)
endif(MSVC)
IF(NOT WIN32)
 IF (REQUIRE_PTHREAD)
  SET(LAPACK_LIBRARY "${LAPACK_LIBRARY};${PTHREAD_LIBRARY}")
  SET(LAPACK_LIBRARIES "${LAPACK_LIBRARIES};${PTHREAD_LIBRARY}")
 ENDIF (REQUIRE_PTHREAD)
ENDIF(NOT WIN32)
if(MAC_VECLIB)
 set(LAPACK_LDFLAGS "-framework vecLib")
endif(MAC_VECLIB) 
set(LAPACK_LINKER_FLAGS ${LAPACK_LDFLAGS})
SET(LAPACK_LIBS ${LAPACK_LIBRARY} ${BLAS_LIBRARY} ${LAPACK_LINKER_FLAGS} CACHE STRING "Flags to link Lapack and Blas (default = ALPS values)")

# FFTW
find_package(FFTW)

# HDF5
# on weiss, it is 2.8.2 and we should not put HL, on 12.04 we need to put it...
if ( ${CMAKE_VERSION} VERSION_LESS "2.8.6") # CHECK THIS BOUND, where are the cmake changelogs ??
 find_package(HDF5 REQUIRED C CXX )
else(${CMAKE_VERSION} VERSION_LESS "2.8.6")
 find_package(HDF5 REQUIRED C CXX HL )
endif(${CMAKE_VERSION} VERSION_LESS "2.8.6")
IF(HDF5_FOUND)
 SET(HAVE_LIBHDF5 1)  
 INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
 triqs_add_definitions(${HDF5_DEFINTIONS})
 SET(ALPS_HAVE_HDF5 1)
ELSE(HDF5_FOUND)
 MESSAGE(FATAL_ERROR "Require hdf5 1.8.2 or higher. Set HDF5_HOME")
ENDIF(HDF5_FOUND)

IF(HDF5_IS_PARALLEL)
 SET(ALPS_HAVE_HDF5_PARALLEL 1)
 MESSAGE(WARNING "parallel(MPI) hdf5 is detected. We will compile but ALPS does not use parallel HDF5. The standard version is preferred.")
 IF(NOT MPI_FOUND)
  MESSAGE(FATAL_ERROR "parallel(MPI) hdf5 needs MPI. Enable MPI or install serial HDF5 libraries.")
 ENDIF(NOT MPI_FOUND)
 INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF(HDF5_IS_PARALLEL)

# python
set(ALPS_BUILD_PYTHON ON)
set(PYTHON_SCRIPTDIR "${CMAKE_INSTALL_PREFIX}/lib/python")
find_package(Python)

IF (PYTHONLIBS_FOUND AND ALPS_BUILD_PYTHON)
 include_directories(SYSTEM ${PYTHON_NUMPY_INCLUDE_DIR})
 MESSAGE (STATUS "Numpy include in ${PYTHON_NUMPY_INCLUDE_DIR}")
 SET(ALPS_HAVE_PYTHON ON)
 INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
 set(BUILD_BOOST_PYTHON TRUE)
ELSE (PYTHONLIBS_FOUND AND ALPS_BUILD_PYTHON)
 set(BUILD_BOOST_PYTHON OFF)
 SET(ALPS_HAVE_PYTHON OFF)
ENDIF (PYTHONLIBS_FOUND AND ALPS_BUILD_PYTHON)

set(ALPS_PYTHON_SITE_PKG ${PYTHON_SITE_PKG})

if (BIND_FORTRAN_LOWERCASE)
 triqs_add_definitions(-DBIND_FORTRAN_LOWERCASE)
endif (BIND_FORTRAN_LOWERCASE)

# Compute various variables

if(MPI_DEFINITIONS)
 set(ALPS_EXTRA_DEFINITIONS "${ALPS_EXTRA_DEFINITIONS} ${MPI_DEFINITIONS}")
endif(MPI_DEFINITIONS)
if(LAPACK_DEFINITIONS)
 set(ALPS_EXTRA_DEFINITIONS "${ALPS_EXTRA_DEFINITIONS} ${LAPACK_DEFINITIONS}")
endif(LAPACK_DEFINITIONS)
if(HDF5_DEFINITIONS)
 set(ALPS_EXTRA_DEFINITIONS "${ALPS_EXTRA_DEFINITIONS} ${HDF5_DEFINITIONS}")
endif(HDF5_DEFINITIONS)

if(MPI_INCLUDE_DIR)
 list(APPEND ALPS_EXTRA_INCLUDE_DIRS ${MPI_INCLUDE_DIR})
 set (MPI_INCLUDE_DIR_FIXED "-I${MPI_INCLUDE_DIR}")
else(MPI_INCLUDE_DIR)
 set (MPI_INCLUDE_DIR_FIXED "")
endif(MPI_INCLUDE_DIR)
#if(SQLite_INCLUDE_DIR)
#  list(APPEND ALPS_EXTRA_INCLUDE_DIRS ${SQLite_INCLUDE_DIR})
#endif(SQLite_INCLUDE_DIR)
if(HDF5_INCLUDE_DIR)
 list(APPEND ALPS_EXTRA_INCLUDE_DIRS ${HDF5_INCLUDE_DIR})
endif(HDF5_INCLUDE_DIR)
if(PYTHON_INCLUDE_DIRS)
 list(APPEND ALPS_EXTRA_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS})
endif(PYTHON_INCLUDE_DIRS)
if(PYTHON_NUMPY_INCLUDE_DIR)
 list(APPEND ALPS_EXTRA_INCLUDE_DIRS ${PYTHON_NUMPY_INCLUDE_DIR})
endif(PYTHON_NUMPY_INCLUDE_DIR)
if(Boost_ROOT_DIR)
 list(APPEND ALPS_EXTRA_INCLUDE_DIRS ${Boost_ROOT_DIR})
endif(Boost_ROOT_DIR)
if(Boost_INCLUDE_DIR)
 list(APPEND ALPS_EXTRA_INCLUDE_DIRS ${Boost_INCLUDE_DIR})
endif(Boost_INCLUDE_DIR)

if(MPI_LINKER_FLAGS)
 set(ALPS_EXTRA_LINKER_FLAGS "${ALPS_EXTRA_LINKER_FLAGS} ${MPI_LINKER_FLAGS}")
endif(MPI_LINKER_FLAGS)
if(LAPACK_LINKER_FLAGS)
 set(ALPS_EXTRA_LINKER_FLAGS "${ALPS_EXTRA_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
endif(LAPACK_LINKER_FLAGS)


if(MPI_LIBRARIES)
 list(APPEND ALPS_EXTRA_LIBRARIES ${MPI_LIBRARIES})
endif(MPI_LIBRARIES)
if(LAPACK_LIBRARIES)
 list(APPEND ALPS_EXTRA_LIBRARIES ${LAPACK_LIBRARIES})
endif(LAPACK_LIBRARIES)
#if(SQLite_LIBRARIES)
#  list(APPEND ALPS_EXTRA_LIBRARIES ${SQLite_LIBRARIES})
#endif(SQLite_LIBRARIES)
#MESSAGE( " HDF5_CXX_LIBRARIES = ${HDF5_CXX_LIBRARIES} ")
#MESSAGE( " HDF5_C_LIBRARIES = ${HDF5_C_LIBRARIES} ")
MESSAGE( STATUS " HDF5_LIBRARIES = ${HDF5_LIBRARIES} ")

if(HDF5_LIBRARIES)
 list(APPEND ALPS_EXTRA_LIBRARIES ${HDF5_LIBRARIES} ) #${HDF5_CXX_LIBRARIES} )
endif(HDF5_LIBRARIES)
if(PYTHON_LIBRARY)
 list(APPEND ALPS_EXTRA_LIBRARIES ${PYTHON_LIBRARY})
endif(PYTHON_LIBRARY)

set (HDF5_LIBS "")
foreach (l ${HDF5_LIBRARIES})
 set (HDF5_LIBS "${HDF5_LIBS} ${l}")
 if (${l} STREQUAL "debug")
  set (HDF5_LIBS "")
 endif()
 if (${l} STREQUAL "optimized")
  set (HDF5_LIBS "")
 endif()
endforeach(l)

if(PYTHON_FOUND)
 set(PYTHON_CPP_FLAGS "-I${PYTHON_INCLUDE_DIRS}")
endif(PYTHON_FOUND)

set (HDF5_LIBRARIES ${HDF5_LIBRARIES_SAVE})

# alps configuration files
set (BIND_FORTRAN_INTEGER_8 ${LAPACK_64_BIT})


# RPATH setting
if(APPLE)
 set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
else(APPLE)
 set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
 set(CMAKE_SKIP_BUILD_RPATH FALSE)
 set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
 set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif(APPLE)

# MPI. 
if (MPI_FOUND)
 triqs_add_definitions( -DHAVE_MPI )
else(MPI_FOUND)
 MESSAGE(FATAL_ERROR "This code requires MPI")
endif(MPI_FOUND)

