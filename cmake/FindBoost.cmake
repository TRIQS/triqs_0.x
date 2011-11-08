

message(STATUS "Boost root dir: ${BOOST_INSTALL_ROOT_DIR}")


include_directories(${BOOST_INSTALL_ROOT_DIR}/include)
set(BOOST_LIBRARY -lboost_python -lboost_serialization -lboost_filesystem -lboost_mpi)

IF(TRIQS_BUILD_STATIC)

  FIND_LIBRARY(BOOST_PYTHON_LIB libboost_python.a ${BOOST_INSTALL_ROOT_DIR}/lib DOC "boost libraries")
  FIND_LIBRARY(BOOST_SERIALIZATION_LIB libboost_serialization.a ${BOOST_INSTALL_ROOT_DIR}/lib DOC "boost libraries")
  FIND_LIBRARY(BOOST_FILESYSTEM_LIB libboost_filesystem.a ${BOOST_INSTALL_ROOT_DIR}/lib DOC "boost libraries")
  FIND_LIBRARY(BOOST_MPI_LIB libboost_mpi.a ${BOOST_INSTALL_ROOT_DIR}/lib DOC "boost libraries")

ELSE(TRIQS_BUILD_STATIC)

  FIND_LIBRARY(BOOST_PYTHON_LIB libboost_python.so ${BOOST_INSTALL_ROOT_DIR}/lib DOC "boost libraries")
  FIND_LIBRARY(BOOST_SERIALIZATION_LIB libboost_serialization.so ${BOOST_INSTALL_ROOT_DIR}/lib DOC "boost libraries")
  FIND_LIBRARY(BOOST_FILESYSTEM_LIB libboost_filesystem.so ${BOOST_INSTALL_ROOT_DIR}/lib DOC "boost libraries")
  FIND_LIBRARY(BOOST_MPI_LIB libboost_mpi.so ${BOOST_INSTALL_ROOT_DIR}/lib DOC "boost libraries")

ENDIF(TRIQS_BUILD_STATIC)
