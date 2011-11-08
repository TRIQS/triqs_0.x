# Compile and link with python 
link_libraries(${PYTHON_LIBRARY})
include_directories(${PYTHON_INCLUDE_DIRS} ${PYTHON_NUMPY_INCLUDE_DIR})

# All the sources
SET (ALLSOURCES 
 ${CMAKE_CURRENT_SOURCE_DIR}/mc_tools/MersenneRNG.cpp
 ${CMAKE_CURRENT_SOURCE_DIR}/mc_tools/polymorphic_random_generator.cpp
 ${CMAKE_CURRENT_SOURCE_DIR}/mc_tools/call_backs.cpp
 ${CMAKE_CURRENT_SOURCE_DIR}/python_tools/python_interpreter.cpp
 ${CMAKE_CURRENT_SOURCE_DIR}/deprecated/utility/my_mpi.cpp 
 ${CMAKE_CURRENT_SOURCE_DIR}/utility/signal_handler.cpp 
 ${CMAKE_CURRENT_SOURCE_DIR}/utility/stack_trace.cpp 
 ${CMAKE_CURRENT_SOURCE_DIR}/deprecated/blitzext/blitz_op.cpp 
)

option(Build_Triqs_General_Tools_Test "Build the test of the general tool ? ( array class , ... ) " ON)
#mark_as_advanced(Build_Array_tests)

# Add the arrays
add_subdirectory(arrays)
add_subdirectory(python_tools)

# The static lib will be build from the sources
add_library(triqs ${ALLSOURCES})
#add_library(triqs SHARED ${ALLSOURCES}) 
target_link_libraries (triqs ${TRIQS_LINK_LIBS} )

# Do we install the dev version?
option(Install_dev "Install the dev version of TRIQS?" OFF)
if (Install_dev)

  # TRIQS cmake file
  configure_file(../cmake/TRIQSConfig.cmake.in Config/TRIQSConfig.cmake)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/Config/TRIQSConfig.cmake DESTINATION share/pytriqs/cmake)

  # Install library
  install(TARGETS triqs DESTINATION lib)

  # Add all headers
  FILE(GLOB ALLHEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    arrays/*.hpp
    arrays/impl/*.hpp
    arrays/IndexMaps/*.hpp
    arrays/IndexMaps/cuboid/*.hpp
    arrays/Storages/*.hpp
    arrays/expressions/*.hpp
    arrays/python/*.hpp
    arrays/hdf5/*.hpp
    arrays/tools/*.hpp
    mc_tools/*.hpp
    python_tools/*.hpp
    utility/*.hpp
  )

  # Install headers in good include dir
  FOREACH(HEADER ${ALLHEADERS})
    get_filename_component(P ${HEADER} PATH)
    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/${HEADER}" DESTINATION "include/triqs/${P}")
  ENDFOREACH(HEADER ${ALLHEADERS})

endif (Install_dev)
