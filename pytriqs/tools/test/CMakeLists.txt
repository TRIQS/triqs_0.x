find_package(TestScriptHdf)
# general tools 
python_build_module(test_converter _pytriqs_test_converter tools/test test_converter.cpp)
add_triqs_test_txt(test_converter )


