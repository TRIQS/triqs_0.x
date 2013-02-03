from gf_imfreq import GfImFreq

cdef class GfImFreq_cython ( GfGeneric_cython ) :
    cdef gf_imfreq _c
    def __init__(self, MeshImFreq mesh, data, TailGf tail, symmetry, indices, name):
        
        GfGeneric_cython.__init__(self, mesh, data, tail, symmetry, indices, name, GfImFreq)
        self._c =  gf_imfreq ( mesh._c, array_view[dcomplex,THREE,COrder](data), tail._c , nothing(), make_c_indices(indices[0],indices[1]) ) 
    
    def __write_hdf5__ (self, gr , char * key) :
        h5_write (make_h5_group_or_file(gr), key, self._c)

    def set_from_fourier(self,GfImTime_cython gt) :
        """Fills self with the Fourier transform of gt"""
        self._c = lazy_fourier( gt._c )

    def set_from_legendre(self, GfLegendre_cython gl) :
        """Fills self with the Legendre transform of gl"""
        self._c = lazy_legendre_imfreq(gl._c)

    def density(self):
        return density(self._c).to_python()

#----------------  Reading from h5 ---------------------------------------

def h5_read_GfImFreq ( gr, std_string key) : 
    return make_GfImFreq( h5_extractor[gf_imfreq]()(make_h5_group_or_file(gr),key))

from pytriqs.base.archive.hdf_archive_schemes import register_class
register_class (GfImFreq, read_fun = h5_read_GfImFreq)

#----------------  Convertions functions ---------------------------------------
  
# Python -> C
cdef gf_imfreq  as_gf_imfreq (g) except +: 
    return (<GfImFreq_cython?>g)._c

# C -> Python. Do NOT add except +
cdef make_GfImFreq ( gf_imfreq  x) :
    return GfImFreq( 
            mesh = make_MeshImFreq (x.mesh()), 
            data = x.data_view().to_python(),
            tail = make_TailGf (x.singularity_view()),
            indices_pack = x.indices()(),
            name = "")

# Python -> C for blocks
cdef gf_block_imfreq  as_gf_block_imfreq (G) except +:
        cdef vector[gf_imfreq] v_c
        for item in G:
            v_c.push_back(as_gf_imfreq(item))
        return make_gf_block_imfreq (v_c)

# C -> Python for block
cdef make_BlockGfImFreq (gf_block_imfreq G) :
    gl = []
    name_list = G.mesh().domain().names()
    for i,n in enumerate(name_list):
        gl.append( make_GfImFreq(G[i] ) )
    return BlockGf( name_list = name_list, block_list = gl)


