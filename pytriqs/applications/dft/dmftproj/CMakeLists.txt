# List the sources
set (SOURCES modules dmftproj readcomline set_ang_trans setsym 
     set_rotloc timeinv read_k_list set_projections orthogonal 
     rot_projectmat density symmetrize_mat rot_dens
     orthogonal_wannier outputqmc outbwin outband)

# add the extension and the path
FOREACH(f ${SOURCES} )
 set(S "${CMAKE_CURRENT_SOURCE_DIR}/${f}.f;${S}")
ENDFOREACH(f)
   
# Name of the program
set(ExecutableName dmftproj)

# The main target and what to link with...
add_executable(${ExecutableName} ${S})
target_link_libraries(${ExecutableName}  ${LAPACK_LIBS} )

# where to install
install (TARGETS ${ExecutableName} DESTINATION bin )

# that is it !


