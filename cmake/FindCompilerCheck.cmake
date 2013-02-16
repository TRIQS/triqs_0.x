#  Copyright Olivier Parcollet 2012 
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

#
#

# first detect the compiler
EXECUTE_PROCESS(COMMAND ${CMAKE_CXX_COMPILER} --version  
 OUTPUT_VARIABLE _compiler_output RESULT_VARIABLE returncode OUTPUT_STRIP_TRAILING_WHITESPACE)

if(CMAKE_COMPILER_IS_GNUCXX )

 SET( compiler_version_min "4.6.3")
 SET( compiler_name "gcc")
 string(REGEX REPLACE ".*([2-5]\\.[0-9]\\.[0-9]).*" "\\1" compiler_version ${_compiler_output})

else(CMAKE_COMPILER_IS_GNUCXX )

 if (_compiler_output MATCHES clang)
  SET(CMAKE_COMPILER_IS_CLANG TRUE )
  SET( compiler_version_min "3.1")
  SET( compiler_name "clang")
  string(REGEX REPLACE ".*([2-5]\\.[0-9]).*" "\\1" compiler_version ${_compiler_output})
 endif (_compiler_output MATCHES clang)

 if (_compiler_output MATCHES Intel)
  SET(CMAKE_COMPILER_IS_ICC TRUE )
  SET( compiler_version_min "13.0.0")
  SET( compiler_name "Intel icc")
  string(REGEX REPLACE "[^0-9]*([0-9]+\\.[0-9]\\.[0-9]).*" "\\1" compiler_version ${_compiler_output})
 endif (_compiler_output MATCHES Intel)

 if (NOT compiler_name) 
  MESSAGE( WARNING "***** Compiler unknown : TRIQS may compile .. or not *******") 
 endif (NOT compiler_name) 

endif(CMAKE_COMPILER_IS_GNUCXX )

MESSAGE( STATUS "Compiler is ${compiler_name} with version ${compiler_version}")

# Check version 
if(compiler_version VERSION_LESS ${compiler_version_min} )
 SET(line_of_star "\n************************** FATAL ERROR ************************************\n")
 MESSAGE( FATAL_ERROR "${line_of_star}You are using the ${compiler_name} compiler but your compiler is too old :\n TRIQS requires version >= ${compiler_version_min} while you have ${compiler_version}\n  ${line_of_star}")
endif(compiler_version VERSION_LESS ${compiler_version_min} )

# Now add some definitions ...

IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
 if (CMAKE_COMPILER_IS_CLANG) 
  add_definitions( -stdlib=libc++ )
 MESSAGE(STATUS " Adding compilation flags -stdlib=libc++ ")
 else (CMAKE_COMPILER_IS_CLANG) 
  MESSAGE( WARNING "${line_of_star}You are on Os X but your are not using clang. This is NOT recommended...${line_of_star}") 
 endif (CMAKE_COMPILER_IS_CLANG) 
ENDIF( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

IF(CMAKE_COMPILER_IS_ICC)
 add_definitions( -DTRIQS_WORKAROUND_INTEL_COMPILER_BUGS)
 add_definitions( -DBOOST_MATH_DISABLE_STD_FPCLASSIFY)
 MESSAGE(STATUS " Adding compilation flags -DTRIQS_WORKAROUND_INTEL_COMPILER_BUGS -DBOOST_MATH_DISABLE_STD_FPCLASSIFY")
ENDIF(CMAKE_COMPILER_IS_ICC)


