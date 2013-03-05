include_directories(${CMAKE_CURRENT_SOURCE_DIR})
cython_module( array_cython array_cython base/test2 )
add_triqs_test_txt(array_cython )

