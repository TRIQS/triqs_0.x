#  Copyright Olivier Parcollet 2013
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

#
# This module looks for triqs
# It sets up : 
#    HERE MAKE THE LIST ! 
# 

SET(TRIAL_PATHS
 $ENV{TRIQS_PATH}/include
 ${TRIQS_PATH}/include
 /usr/include
 /usr/local/include
 /opt/local/include
 /sw/include
 )
FIND_PATH(TRIQS_INCLUDE_DIR triqs_config.h ${TRIAL_PATHS} DOC "Include triqs")

SET(TRIAL_LIBRARY_PATHS
 /usr/lib 
 /usr/local/lib
 /opt/local/lib
 /sw/lib
 $ENV{TRIQS_PATH}/lib
 ${TRIQS_PATH}/lib
 )

SET(TRIQS_LIBRARIES "TRIQS_LIBRARIES-NOTFOUND" CACHE STRING "TRIQS library")
# Try to detect the lib
FIND_LIBRARY(TRIQS_LIBRARIES triqs ${TRIAL_LIBRARY_PATHS} DOC "TRIQS library")

SET(TRIAL_PATHS
 $ENV{TRIQS_PATH}/share/triqs/cmake
 ${TRIQS_PATH}/share/triqs/cmake
 /usr/share/triqs/cmake
 /usr/local/share/triqs/cmake
 /opt/local/share/triqs/cmake
 /sw/share/triqs/cmake
 )
FIND_PATH(TRIQS_CONFIG_FILE TRIQSConfig.cmake ${TRIAL_PATHS} DOC "TRIQS configuration file")

mark_as_advanced(TRIQS_INCLUDE_DIR)
mark_as_advanced(TRIQS_LIBRARIES)

# include prepare variables 
include(${TRIQS_CONFIG_FILE})

FIND_PACKAGE_HANDLE_STANDARD_ARGS(TRIQS DEFAULT_MSG TRIQS_LIBRARIES TRIQS_INCLUDE_DIR)

