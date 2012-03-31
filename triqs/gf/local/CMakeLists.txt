include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)

#add_definitions( -std=gnu++0x) 
add_subdirectory(test)
#add_subdirectory(examples)

function(add_source dir )
foreach (d ${ARGN})
 set_property (GLOBAL APPEND PROPERTY LAZY_SOURCES_LIST ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/${d} )
endforeach (d ${ARGN})
endfunction(add_source dir )

add_source ( . boost_function_lazy.hpp  CMakeLists.txt  expression_container.hpp  lazy_expressions.hpp)
get_property(LAZY_SOURCES_LIST GLOBAL PROPERTY LAZY_SOURCES_LIST)

option (PsPrintCodesArray "Generate a pretty print in ps of the lazy lib" OFF)
mark_as_advanced (PsPrintCodesArray)
if (PsPrintCodesArray)
 find_package(Enscript)
 if (ENSCRIPT_EXECUTABLE)
  SET(code_filename ${CMAKE_BINARY_DIR}/lazy_code.ps)
  add_custom_command (OUTPUT ${code_filename} COMMAND ${ENSCRIPT_EXECUTABLE} -G -f Courier7 --color -Ecpp -o ${code_filename} ${LAZY_SOURCES_LIST} )
  add_custom_target(print_code_lazy ALL DEPENDS ${code_filename})
 endif (ENSCRIPT_EXECUTABLE)

 find_package(Vim)
 if (VIM_EXECUTABLE)
  SET(edit_code_top ${CMAKE_CURRENT_BINARY_DIR}/edit_code.fake)
  add_custom_command (OUTPUT ${edit_code_top} COMMAND ${VIM_EXECUTABLE} -p ${LAZY_SOURCES_LIST} )
  add_custom_target(edit_code_lazy DEPENDS ${edit_code_top})
 endif (VIM_EXECUTABLE)
 
endif (PsPrintCodesArray)


