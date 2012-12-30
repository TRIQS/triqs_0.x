find_package(TestScriptHdf)
add_triqs_test_hdf(wien2k_convert " -p 1.e-6" )
add_triqs_test_hdf(U_mat " -d 1.e-6" )


