# Add C11 flag
SET(CMAKE_CXX_FLAGS "-std=c++0x")

# we do use triqs_array with python support
add_definitions ( -DTRIQS_ARRAYS_WITH_PYTHON_SUPPORT)

install (FILES ${TRIQS_SOURCE_DIR}/pytriqs/__init__.py  DESTINATION ${TRIQS_PYTHON_LIB_DEST} )

add_subdirectory(Tools)
add_subdirectory(Base)
add_subdirectory(Wien2k)
add_subdirectory(Solvers)
add_subdirectory(Solvers/HybridizationExpansion)
add_subdirectory(Solvers/HubbardI)

