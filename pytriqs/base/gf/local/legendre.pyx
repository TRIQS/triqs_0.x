from gf_legendre import GfLegendre

cdef class GfLegendre_cython ( GfGeneric_cython ) :
    cdef gf_legendre _c
    def __init__(self, MeshLegendre mesh, data, tail, symmetry, indices, name):

        GfGeneric_cython.__init__(self, mesh, data, tail, symmetry, indices, name, GfLegendre) 
        self._c =  gf_legendre ( mesh._c, array_view[double,THREE,COrder](data), nothing(), nothing(), make_c_indices(indices[0],indices[1]) ) 
    
    def __write_hdf5__ (self, gr , char * key) :
        h5_write (make_h5_group_or_file(gr), key, self._c)

    def set_from_imtime(self, GfImTime_cython gt) :
        """Fills self with the Legendre transform of gt"""
        self._c = lazy_imtime_legendre(gt._c)

    def set_from_imfreq(self, GfImFreq_cython gw) :
        """Fills self with the Legendre transform of gw"""
        self._c = lazy_imfreq_legendre(gw._c)

    def density(self):
        return density(self._c).to_python()

#----------------  Reading from h5 ---------------------------------------

def h5_read_GfLegendre ( gr, std_string key) : 
    return make_GfLegendre( h5_extractor[gf_legendre]()(make_h5_group_or_file(gr),key))

from pytriqs.base.archive.hdf_archive_schemes import register_class
register_class (GfLegendre, read_fun = h5_read_GfLegendre)

#----------------  Convertions functions ---------------------------------------

# Python -> C
cdef gf_legendre as_gf_legendre (g) except +: 
    return (<GfLegendre_cython?>g)._c

# C -> Python. Do NOT add except +
cdef make_GfLegendre ( gf_legendre  x) :
    return GfLegendre( 
            mesh = make_MeshLegendre (x.mesh()), 
            data = x.data_view().to_python(),
            indices_pack = x.indices()(),
            name = "")

# Python -> C for blocks
cdef gf_block_legendre  as_gf_block_legendre (G) except +:
        cdef vector[gf_legendre] v_c
        for n,g in G:
            v_c.push_back(as_gf_legendre(g))
        return make_gf_block_legendre (v_c)

# C -> Python for block
cdef make_BlockGfLegendre (gf_block_legendre G) :
    gl = []
    name_list = G.mesh().domain().names()
    for i,n in enumerate(name_list):
        gl.append( make_GfLegendre(G[i] ) )
    return BlockGf( name_list = name_list, block_list = gl)

