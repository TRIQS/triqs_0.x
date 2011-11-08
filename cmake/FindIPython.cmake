#
#  Python settings : 
#
#
# The function EXEC_PYTHON_SCRIPT executes the_script in  python interpreter
# and set the variable of output_var_name in the calling scope
#
FUNCTION ( EXEC_PYTHON_SCRIPT the_script output_var_name)
  EXECUTE_PROCESS(COMMAND ${PYTHON_INTERPRETER} -c "${the_script}" 
    OUTPUT_VARIABLE res RESULT_VARIABLE returncode OUTPUT_STRIP_TRAILING_WHITESPACE)
#  IF (NOT returncode EQUAL 0)
#    MESSAGE(FATAL_ERROR "The script : ${the_script} \n did not run properly in the Python interpreter. Check your python installation.") 
#  ENDIF (NOT returncode EQUAL 0)
  SET( ${output_var_name} ${res} PARENT_SCOPE)
ENDFUNCTION (EXEC_PYTHON_SCRIPT)

EXEC_PYTHON_SCRIPT ("import distutils " nulle) # check that distutils is there...
EXEC_PYTHON_SCRIPT ("import numpy" nulle) # check that numpy is there...
#EXEC_PYTHON_SCRIPT ("import tables" nulle) # check that tables is there...
#EXEC_PYTHON_SCRIPT ("import h5py" nulle) # check that h5py is there...
#EXEC_PYTHON_SCRIPT ("import pyalps.hdf5" nulle) # check that pyalps.hdf5 is there...
MESSAGE(STATUS "Python interpreter ok")

#EXEC_PYTHON_SCRIPT ("try :\n import IPython.Shell\n print 1\nexcept:\n print 0" IPYTHON_AVAIL) # check that ipython is there...
EXEC_PYTHON_SCRIPT ("try :\n import IPython\n print 1\nexcept:\n print 0" IPYTHON_AVAIL) # check that ipython is there...


