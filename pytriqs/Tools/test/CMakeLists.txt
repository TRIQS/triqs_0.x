find_package(TestScriptHdf)
# general tools 
add_definitions(-DTRIQS_ARRAYS_WITH_PYTHON_SUPPORT)
python_build_module(test_converter _pytriqs_test_converter Tools/test test_converter.cpp)
add_triqs_test_txt(test_converter )


