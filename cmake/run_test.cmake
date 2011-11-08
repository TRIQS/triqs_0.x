set (output_new  ${CMAKE_CURRENT_BINARY_DIR}/${name}_output) 

# hum, only a  recent version of hdf5 1.8.7 ?
#SET(H5_DIFF_EXCLUDE --exclude-path=revisions  --exclude-path=log_type --exclude-path=state_type)

if (H5_DIFF_EXECUTABLE)
 separate_arguments(H5_DIFF_OPTIONS)
 separate_arguments(H5_DIFF_OBJECTS)
 set (COM  ${H5_DIFF_EXECUTABLE} ${H5_DIFF_EXCLUDE}  ${H5_DIFF_OPTIONS} ${outputName}  ${reference} ${H5_DIFF_OBJECTS})
else(H5_DIFF_EXECUTABLE)
 set (COM  ${CMAKE_COMMAND} -E compare_files ${output_new} ${reference})  
endif (H5_DIFF_EXECUTABLE)

message(" command for the test ${cmd} ${input}")
execute_process(
 COMMAND ${cmd}
 RESULT_VARIABLE not_successful
 INPUT_FILE ${input}
 OUTPUT_FILE ${output_new}
 ERROR_FILE ${output_new}.err
 ERROR_VARIABLE err
 TIMEOUT 600
 )

if(not_successful)
 message(SEND_ERROR "error runing test '${name}': ${err};shell output: ${not_successful}!")
endif(not_successful)

MESSAGE( "ABOUT TO compare with ${COM}")

execute_process(
 COMMAND ${COM}
 RESULT_VARIABLE not_successful
 OUTPUT_VARIABLE out
 ERROR_VARIABLE err
 TIMEOUT 600
 )

if(not_successful)
 message(SEND_ERROR "output does not match for '${name}': ${err}; ${out}; shell output: ${not_successful}!")
endif(not_successful)
endif(output)

#file(REMOVE ${output_new})

