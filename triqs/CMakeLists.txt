# Add C11 flag
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

# Compile and link with python 
link_libraries(${PYTHON_LIBRARY})
include_directories(${PYTHON_INCLUDE_DIRS} ${PYTHON_NUMPY_INCLUDE_DIR})

option(Use_Pade "Compile Pade tools?" OFF)

# All the sources
SET (ALLSOURCES 
 ${CMAKE_CURRENT_SOURCE_DIR}/mc_tools/MersenneRNG.cpp
 ${CMAKE_CURRENT_SOURCE_DIR}/mc_tools/random_generator.cpp
 ${CMAKE_CURRENT_SOURCE_DIR}/python_tools/python_interpreter.cpp
 ${CMAKE_CURRENT_SOURCE_DIR}/deprecated/utility/my_mpi.cpp 
 ${CMAKE_CURRENT_SOURCE_DIR}/utility/signal_handler.cpp 
 ${CMAKE_CURRENT_SOURCE_DIR}/utility/callbacks.cpp
 ${CMAKE_CURRENT_SOURCE_DIR}/utility/stack_trace.cpp 
 ${CMAKE_CURRENT_SOURCE_DIR}/deprecated/blitzext/blitz_op.cpp 
 ${CMAKE_CURRENT_SOURCE_DIR}/gf_local/fourier.cpp
 ${CMAKE_CURRENT_SOURCE_DIR}/gf_local/GF_Bloc_Base.cpp
 ${CMAKE_CURRENT_SOURCE_DIR}/gf_local/GF_Bloc_ImFreq.cpp
 ${CMAKE_CURRENT_SOURCE_DIR}/gf_local/GF_Bloc_ImLegendre.cpp
 ${CMAKE_CURRENT_SOURCE_DIR}/gf_local/GF_Bloc_ImTime.cpp
 ${CMAKE_CURRENT_SOURCE_DIR}/gf_local/GF_Bloc_ReFreq.cpp
 ${CMAKE_CURRENT_SOURCE_DIR}/gf_local/GF_Bloc_ReTime.cpp
 ${CMAKE_CURRENT_SOURCE_DIR}/gf_local/legendre_matsubara.cpp
 ${CMAKE_CURRENT_SOURCE_DIR}/gf_local/MeshGF.cpp
 ${CMAKE_CURRENT_SOURCE_DIR}/gf_local/TailGF.cpp

 ${CMAKE_CURRENT_SOURCE_DIR}/gf/local/fourier_base.cpp
 ${CMAKE_CURRENT_SOURCE_DIR}/gf/local/fourier_matsubara.cpp
 # ${CMAKE_CURRENT_SOURCE_DIR}/gf/local/fourier_real.cpp
 ${CMAKE_CURRENT_SOURCE_DIR}/gf/local/density.cpp

)

if (Use_Pade)
  set (ALLSOURCES ${ALLSOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/gf_local/pade.cpp)
  add_definitions(-DUSE_PADE)
  find_package(GMP)
  SET (TRIQS_LINK_LIBS ${TRIQS_LINK_LIBS} ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})
endif (Use_Pade)

option(Build_Triqs_General_Tools_Test "Build the test of the general tool ? ( array class , ... ) " ON)

# Add the arrays
add_subdirectory(utility)
add_subdirectory(arrays)
add_subdirectory(python_tools)
add_subdirectory(clef)
add_subdirectory(det_manip)
add_subdirectory(gf)

# The lib will be build from the sources
add_library(triqs ${ALLSOURCES})
target_link_libraries (triqs ${TRIQS_LINK_LIBS} )

# Install library
install(TARGETS triqs DESTINATION lib)

# Do we install the dev version?
option(Install_dev "Install the headers and libraries of TRIQS?" ON)
if (Install_dev)

  # TRIQS cmake file
  configure_file(../cmake/TRIQSConfig.cmake.in Config/TRIQSConfig.cmake)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/Config/TRIQSConfig.cmake DESTINATION share/triqs/cmake)


  # Add all headers
  FILE(GLOB ALLHEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    arrays/*.hpp
    arrays/storages/*.hpp
    arrays/python/*.hpp
    arrays/impl/*.hpp
    arrays/functional/*.hpp
    arrays/proto/*.hpp
    arrays/h5/*.hpp
    arrays/linalg/*.hpp
    arrays/indexmaps/*.hpp
    arrays/indexmaps/cuboid/*.hpp
    mc_tools/*.hpp
    python_tools/*.hpp
    gf_local/*.hpp
    utility/*.hpp
    utility/proto/*.hpp
    clef/*.hpp
    clef/adapters/*.hpp
  )

  # Install headers in good include dir
  FOREACH(HEADER ${ALLHEADERS})
    get_filename_component(P ${HEADER} PATH)
    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/${HEADER}" DESTINATION "include/triqs/${P}")
  ENDFOREACH(HEADER ${ALLHEADERS})

endif (Install_dev)
