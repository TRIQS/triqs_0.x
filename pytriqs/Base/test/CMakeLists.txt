find_package(TestScriptHdf)
if (MPI_FOUND)
 find_package(MPI)
 MESSAGE(STATUS "MPI launch command is : ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} PROCS ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
endif (MPI_FOUND)

# A sample test with text output
add_triqs_test_txt(ExampleTest )

# A sample test with hdf5 output
add_triqs_test_hdf(ExampleTestH5 " -p 1.e-6" )
#add_triqs_test_hdf(ExampleTestH5-2 " -p 1.e-6" )

# Some basic HDF IO stuff:
add_triqs_test_hdf(HDF5_IO " -p 1.e-6" )

# Greens functions:
add_triqs_test_hdf(GF_Init " -p 1.e-6" )
add_triqs_test_hdf(GF_BasOps " -p 1.e-6" ) #" -p 1.e-6" ) # quick fix for pb on Ubuntu... big relative diff ?

# a simple dos on square lattice
add_triqs_test_hdf(dos1 " -d 1.e-6" )

# Pade approximation
add_triqs_test_hdf(Pade " -d 1.e-6" )

