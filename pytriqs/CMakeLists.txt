# Create a temporary copy of the python modules so that we can run before installation with the test
#FILE(COPY ${TRIQS_SOURCE_DIR}/pytriqs DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ FILES_MATCHING PATTERN *.py PATTERN *.svn EXCLUDE )
FILE(GLOB_RECURSE all_py_files RELATIVE ${CMAKE_SOURCE_DIR}/pytriqs *.py )
SET(py_copy_tar  ${CMAKE_BINARY_DIR}/py_copy.tar)
add_custom_command (OUTPUT ${py_copy_tar} DEPENDS ${all_py_files}
 COMMAND cd ${CMAKE_SOURCE_DIR}/pytriqs && tar cf ${py_copy_tar} ${all_py_files} && cd ${CMAKE_BINARY_DIR}/pytriqs && tar xf ${py_copy_tar} )
add_custom_target(py_copy ALL DEPENDS ${py_copy_tar})

# 
find_package(CythonModuleMacro)

# Add C11 flag
# Put this when boost python will be removed (or fixed for C++11).
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

# for CYTHON ONLY
add_definitions (-DTRIQS_ALLOW_EMPTY_VIEW)

install (FILES ${TRIQS_SOURCE_DIR}/pytriqs/__init__.py  DESTINATION ${TRIQS_PYTHON_LIB_DEST} )

add_all_subdirectories_with_cmakelist()

