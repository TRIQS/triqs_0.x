
find_package(TestCSimple)
find_package(TestScriptHdf)

# Build python module
SET(SOURCES_WRAP test1.cpp)
python_build_module(ArrayTests  _array_tests base/Tests ${SOURCES_WRAP})

if (NOT TRIQS_BUILD_STATIC)
 target_link_libraries(_array_tests  triqs ${TRIQS_LINK_LIBS})
endif (NOT TRIQS_BUILD_STATIC)


# A sample test with text output
add_triqs_test_txt( test1 )

