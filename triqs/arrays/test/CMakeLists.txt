set (Build_Array_tests ${Build_Triqs_General_Tools_Test})
if (Build_Array_tests)
add_subdirectory (C++)
add_subdirectory (C++_with_py_support)
#add_subdirectory (python)
#add_subdirectory (speed)
endif (Build_Array_tests)

