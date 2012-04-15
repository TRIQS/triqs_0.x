include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)

add_subdirectory(test)

function(add_source dir )
foreach (d ${ARGN})
 set_property (GLOBAL APPEND PROPERTY ARRAY_SOURCES_LIST ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/${d} )
endforeach (d ${ARGN})
endfunction(add_source dir )

add_source ( . array.hpp vector.hpp matrix.hpp)
add_source ( impl option.hpp sliceable_object.hpp indexmap_storage_pair.hpp assignment.hpp compound_assign.hpp iterator_adapter.hpp cache.hpp common.hpp tuple_tools.hpp )
add_source ( storages mem_block.hpp shared_block.hpp )
add_source ( python numpy_interface.hpp  converters.hpp)
add_source ( h5 common.hpp simple_read_write.hpp array_stack.hpp array_proxy.hpp )
add_source ( indexmaps/cuboid cuboid_domain.hpp index_order.hpp cuboid_map.hpp cuboid_index_generator.hpp cuboid_map_iterator.hpp foreach.hpp cuboid_slice.hpp /../permutation.hpp )
add_source ( expressions arithmetic.hpp algebra_common.hpp array_algebra.hpp matrix_algebra.hpp map.hpp fold.hpp  )
add_source ( linalg matrix_cache.hpp a_x_ty.hpp det_and_inverse.hpp matmul.hpp mat_vec_mul.hpp eigenelements.hpp )
get_property(ARRAY_SOURCES_LIST GLOBAL PROPERTY ARRAY_SOURCES_LIST)

option (PsPrintCodesArray "Generate a pretty print in ps of the array lib" OFF)
mark_as_advanced (PsPrintCodesArray)
if (PsPrintCodesArray)
 find_package(Enscript)
 if (ENSCRIPT_EXECUTABLE)
  SET(code_filename ${CMAKE_BINARY_DIR}/array_code.ps)
  add_custom_command (OUTPUT ${code_filename} COMMAND ${ENSCRIPT_EXECUTABLE} -G -f Courier7 --color -Ecpp -o ${code_filename} ${ARRAY_SOURCES_LIST} )
  add_custom_target(print_code_array ALL DEPENDS ${code_filename})
 endif (ENSCRIPT_EXECUTABLE)

 find_package(Vim)
 if (VIM_EXECUTABLE)
  SET(edit_code_top ${CMAKE_CURRENT_BINARY_DIR}/edit_code.fake)
  add_custom_command (OUTPUT ${edit_code_top} COMMAND ${VIM_EXECUTABLE} -p ${ARRAY_SOURCES_LIST} )
  add_custom_target(edit_code_array DEPENDS ${edit_code_top})
 endif (VIM_EXECUTABLE)

endif (PsPrintCodesArray)


