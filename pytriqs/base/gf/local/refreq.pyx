from gf_refreq import GfReFreq

cdef class GfReFreq_cython ( GfGeneric_cython ) :
    cdef gf_refreq _c
    def __init__(self, MeshReFreq mesh, data, TailGf tail, symmetry, indices,  name ):

        GfGeneric_cython.__init__(self,  mesh, data,  tail, symmetry,indices, name, GfReFreq) 
        self._c =  gf_refreq ( mesh._c, array_view[dcomplex,THREE](data), tail._c , nothing(), make_c_indices(indices[0],indices[1]) ) 
    
    def __write_hdf5__ (self, gr , char * key) :
        h5_write (make_h5_group(gr), key, self._c)

    def set_from_pade(self, GfImFreq_cython gw, n_matsubara_freq = 100, freq_offset = 0.0) :
        pade(self._c, gw._c, n_matsubara_freq, freq_offset)

#----------------  Reading from h5 ---------------------------------------

def h5_read_GfReFreq ( gr, std_string key) : 
    return make_GfReFreq( h5_extractor[gf_refreq]()(make_h5_group(gr),key))

from pytriqs.base.archive.hdf_archive_schemes import register_class
register_class (GfReFreq, read_fun = h5_read_GfReFreq)

#----------------  Convertions functions ---------------------------------------

# Python -> C
cdef gf_refreq as_gf_refreq (g) except +:
    return (<GfReFreq_cython?>g)._c

# C -> Python. Do NOT add except +
cdef make_GfReFreq ( gf_refreq  x) :
    return GfReFreq( 
            mesh = make_MeshReFreq (x.mesh()), 
            data = x.data_view().to_python(),
            tail = make_TailGf (x.singularity_view()),
            indices_pack = x.indices()(),
            name = "")

# Python -> C for blocks
cdef gf_block_refreq  as_gf_block_refreq (G) except +:
        cdef vector[gf_refreq] v_c
        for n,g in G:
            v_c.push_back(as_gf_refreq(g))
        return make_gf_block_refreq (v_c)

# C -> Python for block
cdef make_BlockGfReFreq (gf_block_refreq G) :
    gl = []
    name_list = G.mesh().domain().names()
    for i,n in enumerate(name_list):
        gl.append( make_GfReFreq(G[i] ) )
    return BlockGf( name_list = name_list, block_list = gl)

