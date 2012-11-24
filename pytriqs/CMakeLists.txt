# 
find_package(CythonModuleMacro)

# Add C11 flag
# Put this when boost python will be removed (or fixed for C++11).
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

# for CYTHON ONLY
add_definitions (-DTRIQS_ALLOW_EMPTY_VIEW)

install (FILES ${TRIQS_SOURCE_DIR}/pytriqs/__init__.py  DESTINATION ${TRIQS_PYTHON_LIB_DEST} )

add_all_subdirectories_with_cmakelist()

