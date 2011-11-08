#  Copyright Olivier Parcollet 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)
#

FIND_PROGRAM(VIM_EXECUTABLE
 NAMES gvim mvim vim vi  
 PATHS /usr/bin /opt/local/bin $ENV{HOME}/bin
 PATH_SUFFIXES  bin
 )

if (NOT VIM_EXECUTABLE)
 MESSAGE(STATUS "Info: I can not find gvim, mvim, vim or vi: the edit_code targets will not work [developers only]")
else (NOT VIM_EXECUTABLE)
 INCLUDE(FindPackageHandleStandardArgs)
 FIND_PACKAGE_HANDLE_STANDARD_ARGS(EDIT_CODE DEFAULT_MSG VIM_EXECUTABLE)
 MARK_AS_ADVANCED( VIM_EXECUTABLE )
endif (NOT VIM_EXECUTABLE)


