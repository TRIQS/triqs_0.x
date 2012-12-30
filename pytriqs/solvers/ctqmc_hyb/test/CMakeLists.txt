find_package(TestScriptHdf)

add_triqs_test_hdf(SingleSiteBethe " -p 1.e-5" )
add_triqs_test_hdf(CDMFT_4_sites " -p 1.e-5"  )


