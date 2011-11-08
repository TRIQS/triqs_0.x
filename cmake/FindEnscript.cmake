#  Copyright Olivier Parcollet 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)
#

FIND_PROGRAM(ENSCRIPT_EXECUTABLE
 NAMES enscript 
 PATHS /usr/bin /opt/local/bin $ENV{HOME}/bin
 PATH_SUFFIXES  bin
 )

if (NOT ENSCRIPT_EXECUTABLE)
 MESSAGE(STATUS "Info: I can not enscript: the print_code targets will not work [developers only]")
else (NOT ENSCRIPT_EXECUTABLE)
 INCLUDE(FindPackageHandleStandardArgs)
 FIND_PACKAGE_HANDLE_STANDARD_ARGS(PRINT_CODE DEFAULT_MSG ENSCRIPT_EXECUTABLE)
 MARK_AS_ADVANCED( ENSCRIPT_EXECUTABLE )
endif (NOT ENSCRIPT_EXECUTABLE)



