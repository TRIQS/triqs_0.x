#  Copyright Olivier Parcollet 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

# This forbids compilation in the source directory

IF (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  MESSAGE(FATAL_ERROR "I am sorry, the on-site compilation is disabled at the moment. Use another directory")
ENDIF (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})


