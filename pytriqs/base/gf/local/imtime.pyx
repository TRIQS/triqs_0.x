from gf_imtime import GfImTime

cdef class GfImTime_cython ( GfGeneric_cython ) :
    cdef gf_imtime _c
    def __init__(self, MeshImTime mesh, data, TailGf tail, symmetry, indices,  name ):

        GfGeneric_cython.__init__(self,  mesh, data,  tail, symmetry,indices, name, GfImTime) 
        self._c =  gf_imtime ( mesh._c, array_view[double,THREE,COrder](data), tail._c , nothing(), make_c_indices(indices[0],indices[1]) ) 
    
    def __write_hdf5__ (self, gr , char * key) :
        h5_write (make_h5_group_or_file(gr), key, self._c)

    def set_from_inverse_fourier_of(self,GfImFreq_cython gw) :
        """Fills self with the Inverse Fourier transform of gw"""        
        self._c = lazy_inverse_fourier( gw._c)

#----------------  Reading from h5 ---------------------------------------

def h5_read_GfImTime ( gr, std_string key) : 
    return make_GfImTime( h5_extractor[gf_imtime]()(make_h5_group_or_file(gr),key))

from pytriqs.base.archive.hdf_archive_schemes import register_class
register_class (GfImTime, read_fun = h5_read_GfImTime)

#----------------  Convertions functions ---------------------------------------

# Python -> C
cdef gf_imtime  as_gf_imtime (g) except +: 
    return (<GfImTime_cython?>g)._c

# C -> Python. Do NOT add except +
cdef make_GfImTime ( gf_imtime  x) :
    return GfImTime( 
            Mesh = make_MeshImTime (x.mesh()), 
            Data = x.data_view().to_python(),
            Tail = make_TailGf (x.singularity_view()),
            IndicesPack = x.indices()(),
            Name = "")

# Python -> C for blocks
cdef gf_block_imtime  as_gf_block_imtime (G) except +:
        cdef vector[gf_imtime] v_c
        for item in G:
            v_c.push_back(as_gf_imtime(item))
        return make_gf_block_imtime (v_c)

# C -> Python for block
cdef make_BlockGfImTime (gf_block_imtime G) :
    gl = []
    name_list = G.mesh().domain().names()
    for i,n in enumerate(name_list):
        gl.append( make_GfImTime(G[i] ) )
    return BlockGf( name_list = name_list, block_list = gl)

