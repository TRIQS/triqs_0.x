find_package(TestScriptHdf)
# general tools 
python_build_module(boost_converter _pytriqs_boost_converter tools/test boost_converter.cpp)
add_triqs_test_txt(boost_converter )

cython_module( array_cython array_cython base/test2 )
add_triqs_test_txt(array_cython )

