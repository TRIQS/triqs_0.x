find_package(TestScriptHdf)
FILE(COPY SrVO3.h5 SrVO3.ctqmcout SrVO3.symqmc SrVO3.sympar SrVO3.parproj DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_triqs_test_hdf(wien2k_convert " -p 1.e-6" )
add_triqs_test_hdf(sumklda_basic " -p 1.e-6" )
add_triqs_test_hdf(U_mat " -d 1.e-6" )


