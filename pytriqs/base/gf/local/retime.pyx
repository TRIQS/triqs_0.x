from gf_retime import GfReTime

cdef class GfReTime_cython ( GfGeneric_cython ) :
    cdef gf_retime _c
    def __init__(self, MeshReTime mesh, data, TailGf tail, symmetry, indices,  name ):

        GfGeneric_cython.__init__(self,  mesh, data,  tail, symmetry,indices, name, GfReTime) 
        self._c =  gf_retime ( mesh._c, array_view[dcomplex,THREE](data), tail._c , nothing(), make_c_indices(indices[0],indices[1]) ) 
    
    def __write_hdf5__ (self, gr , char * key) :
        h5_write (make_h5_group(gr), key, self._c)

#----------------  Reading from h5 ---------------------------------------

def h5_read_GfReTime ( gr, std_string key) : 
    return make_GfReTime( h5_extractor[gf_retime]()(make_h5_group(gr),key))

from pytriqs.base.archive.hdf_archive_schemes import register_class
register_class (GfReTime, read_fun = h5_read_GfReTime)

#----------------  Convertions functions ---------------------------------------

# Python -> C
cdef gf_retime as_gf_retime (g) except +:
    return (<GfReTime_cython?>g)._c

# C -> Python. Do NOT add except +
cdef make_GfReTime ( gf_retime  x) :
    return GfReTime( 
            mesh = make_MeshReTime (x.mesh()), 
            data = x.data_view().to_python(),
            tail = make_TailGf (x.singularity_view()),
            indices_pack = x.indices()(),
            name = "")

# Python -> C for blocks
cdef gf_block_retime  as_gf_block_retime (G) except +:
        cdef vector[gf_retime] v_c
        for n,g in G:
            v_c.push_back(as_gf_retime(g))
        return make_gf_block_retime (v_c)

# C -> Python for block
cdef make_BlockGfReTime (gf_block_retime G) :
    gl = []
    name_list = G.mesh().domain().names()
    for i,n in enumerate(name_list):
        gl.append( make_GfReTime(G[i] ) )
    return BlockGf( name_list = name_list, block_list = gl)

