#  Copyright Olivier Parcollet 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)
#

option (PsPrintCodes "Generate a pretty print in ps of the codes" OFF)
mark_as_advanced (PsPrintCodes)

if (PsPrintCodes)

 if (NOT ENSCRIPT_EXECUTABLE)
  find_package(Enscript)
 endif (NOT ENSCRIPT_EXECUTABLE)

 if (ENSCRIPT_EXECUTABLE)
  SET(code_filenameh ${CMAKE_BINARY_DIR}/${CODENAME}_code_h.ps)
  SET(code_filenamec ${CMAKE_BINARY_DIR}/${CODENAME}_code_c.ps)
  add_custom_command (OUTPUT ${code_filenameh} COMMAND ${ENSCRIPT_EXECUTABLE} -G -f Courier7 --color -Ecpp -o ${code_filenameh} ${SOURCES_HPP} )
  add_custom_command (OUTPUT ${code_filenamec} COMMAND ${ENSCRIPT_EXECUTABLE} -G -f Courier7 --color -Ecpp -o ${code_filenamec} ${SOURCES_CPP} )
  add_custom_target(print_code_${CODENAME} ALL DEPENDS ${code_filenameh} ${code_filenamec} )
 endif (ENSCRIPT_EXECUTABLE)

 if (NOT VIM_EXECUTABLE)
  find_package(Vim)
 endif (NOT VIM_EXECUTABLE)

 if (VIM_EXECUTABLE)
  SET(edit_code_top ${CMAKE_CURRENT_BINARY_DIR}/edit_code.fake)
  add_custom_command (OUTPUT ${edit_code_top} COMMAND ${VIM_EXECUTABLE} -p ${SOURCES_CPP} ${SOURCES_HPP} )
  add_custom_target(edit_code_${CODENAME} DEPENDS ${edit_code_top})
 endif (VIM_EXECUTABLE)

endif (PsPrintCodes)

