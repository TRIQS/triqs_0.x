# 
find_package(CythonModuleMacro)

# Add C11 flag
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

# for CYTHON ONLY
add_definitions (-DTRIQS_ALLOW_EMPTY_VIEW)

install (FILES ${TRIQS_SOURCE_DIR}/pytriqs/__init__.py  DESTINATION ${TRIQS_PYTHON_LIB_DEST} )

add_subdirectory(Tools)
add_subdirectory(Base)
add_subdirectory(Wien2k)
add_subdirectory(Solvers)

add_subdirectory(demo)

