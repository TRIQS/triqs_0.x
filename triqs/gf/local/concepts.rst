
.. highlight:: c

Here are the concepts for a function of one variable.

Concepts
=============================================================

ResultOf 
---------------------

It is the boost::result_of protocol, mainly useful for pre-C++11.

* **Definition** : 

+---------------------------------------------------------------+---------------------------------------------------------------------+
| Elements                                                      | Comment                                                             |
+===============================================================+=====================================================================+
| * template<typename T> struct result                          |                                                                     |
+---------------------------------------------------------------+---------------------------------------------------------------------+
| * template<typename This, typename X> struct result<This(X)>; | Metafunction to compute result type of operator()(X const &x) const |
+---------------------------------------------------------------+---------------------------------------------------------------------+
 

LazyCallable
---------------------

The () operators of the objects can also be called with triqs::lazy expressions.


Domain
------------------------------------------------- 

* **Purpose**  : The domain of definition of a function. It is a mathematical definition of the domain,
  and does not contain any mesh.

* **Refines** : CopyConstructible, DefaultContructible, EqualComparable, BoostSerializable
   ?? Equal ??

* **Definition** : 

+----------------------------------------------------------------------------+---------------------------------------------------------------------+
| Elements                                                                   | Comment                                                             |
+============================================================================+=====================================================================+
| * point_type                                                               | Type of element in the domain (int, int, double, k_vector, ...) as  |
|                                                                            | in the call of function. In particular, in Matsubara, it is int.    |
+----------------------------------------------------------------------------+---------------------------------------------------------------------+
| * embedded_point_type                                                      | Type of element in the domain (int, complex, double, k_vector, ...) |
|                                                                            | The real type in the domain. In particular, in Matsubara, it is a   |
|                                                                            | complex.                                                            |
+----------------------------------------------------------------------------+---------------------------------------------------------------------+

The
* **Examples** :
  
   * Matsubara frequencies (boson/fermion)
   * Matsubara time
   * Real frequencies
   * Real time 
   * Brillouin zone

Mesh
------------------------------------------------- 

* **Purpose**  : A mesh over a domain, and more generally the practical representation of the domain.
  It does not really need to be a mesh : e.g. if the function is represented on a polynomial basis, 
  it is the parameters of this representation (max number of coordinates, e.g.)

* **Refines** : CopyConstructible, DefaultContructible, EqualComparable

* **Examples** : Some domains and the corresponding possible meshes.

+-----------------------------------------------------+--------------------------------------------------------+
| Domain Type                                         | Possible mesh types : what does it contain ?           |
+=====================================================+========================================================+
| * Matsubara frequencies                             | Nmax, the max number of Matsubara Freq.                |
+-----------------------------------------------------+--------------------------------------------------------+
| * Matsubara time                                    | The time slicing is on a mesh, or the Legendre mesh is |
|                                                     | we store the function with Legendre representation.    |
+-----------------------------------------------------+--------------------------------------------------------+
| * Real frequencies                                  | Parameters of the mesh in frequency                    |
+-----------------------------------------------------+--------------------------------------------------------+
| * Real time                                         |                                                        |
+-----------------------------------------------------+--------------------------------------------------------+
| * Brillouin zone                                    | parameters of the mesh over the BZ, symmetry ?         |
+-----------------------------------------------------+--------------------------------------------------------+


* **Definition** : 

+-----------------------------------------------------------------------+------------------------------------------------------------+
| Elements                                                              | Comment                                                    |
+=======================================================================+============================================================+
| * domain_type                                                         | Type of the Domain represented                             |
+-----------------------------------------------------------------------+------------------------------------------------------------+
| * domain_type const & domain() const                                  | Returns the domain                                         |
+-----------------------------------------------------------------------+------------------------------------------------------------+
| * index_type                                                          | Type of indices. e.g. triqs::arrays::range                 |
+-----------------------------------------------------------------------+------------------------------------------------------------+
| * mesh_pt<THIS> operator[](index_type) const                          | From an index, return a mesh_pt (Cf below) containing this |
|                                                                       | a ref to this mesh and the index.                          |
+-----------------------------------------------------------------------+------------------------------------------------------------+
| * domain_type::point_type embed(index_type const &) const           | From the index, return the corresponding point_type      |
+-----------------------------------------------------------------------+------------------------------------------------------------+
| * size_t size() const                                                 | The length of the array needed to store the mesh           |
+-----------------------------------------------------------------------+------------------------------------------------------------+
| * interpolate(F,domain_type::element_type p,Tag=Default)              | Interpolation of the function F on the point p. Tag is     |
|                                                                       | encoding the method to use (with a default).               |
+-----------------------------------------------------------------------+------------------------------------------------------------+
| * slice_args_type                                                     | Type of the argument of slice method                       |
+-----------------------------------------------------------------------+------------------------------------------------------------+
| * Mesh slice(slice_args_type) const                                   | slice method of the mesh                                   |
+-----------------------------------------------------------------------+------------------------------------------------------------+

mesh_pt 
^^^^^^^^^^^^^^^^

mesh_pt is a little wrapper that contains basically a tuple of a reference to a mesh, and the corresponding indice ::

  template<typename MeshType> struct mesh_pt { 
   MeshType const & m; 
   typename MeshType::index_type i;
   mesh_pt( MeshType const & mesh, typename MeshType::index_type const & index): m(mesh), i(index) {}
   
   typedef typename MeshType::domain_type::element_type cast_type;
   operator cast_type() const; // cast into the element type of the domain (e.g. real time, real frequency).

  };

mesh_range
^^^^^^^^^^^^^^^^

mesh_range is a little wrapper that contains basically a tuple of a reference to a mesh, and the corresponding indice ::

  template<typename MeshType> struct mesh_range { 
   MeshType const & m; 
   typename MeshType::range_type r;
   mesh_pt( MeshType const & mesh, typename MeshType::range_typeindex_type const & index): m(mesh), i(index) {}
   
   typedef ??? cast_type;      // cast into a range of element_type : the precise type depends on the  
   operator cast_type() const; // mesh type.

  };

.. note:: Should I unify these class ?? --> not really, very simple . It is clearer like this.


PureFunctionOnDomain 
-----------------------

* **Purpose**  : A function from a domain to a target space. 

* **Refines**   :

* **Definition** : 

+----------------------------------------------+---------------------------------------------------------+
| Elements                                     | Comment                                                 |
+==============================================+=========================================================+
| domain_type                                  | Type of the Domain represented                          |
+----------------------------------------------+---------------------------------------------------------+
| domain_type const & domain() const           | Returns the domain                                      |
+----------------------------------------------+---------------------------------------------------------+
| shape_type                                   | result of shape (mini_vector<size_t,2>)                 |
+----------------------------------------------+---------------------------------------------------------+
| shape_type shape() const                     | Shape of the tail                                       |
+----------------------------------------------+---------------------------------------------------------+
| operator (domain_type::element_type) const   | Calling for all elements of the Domain (including infty |
|                                              | if it is in the domain...                               |
+----------------------------------------------+---------------------------------------------------------+

PureFunctionOnMesh 
-----------------------

* **Purpose**  : A function from a domain to a target space, represented in a mesh. 

* **Refines**   : PureFunctionOnDomain.

* **Definition** : 

+------------------------------------------------+-------------------------------------------------------+
| Elements                                       | Comment                                               |
+================================================+=======================================================+
| * mesh_type                                    | Type of the mesh representing the domain.             |
+------------------------------------------------+-------------------------------------------------------+
| * mesh_type const & mesh() const               | Returns the mesh.                                     |
+------------------------------------------------+-------------------------------------------------------+
| * operator ( grid_pt<mesh_type> const &) const | Calling on a grid_pt gives direct access to the value |
|                                                | on a grid point.                                      |
+------------------------------------------------+-------------------------------------------------------+

NB : the result type of the () operator are either deduces by modeling ResultOf or using C++11 technique, simply...


LocalGf : the immutable local Green function
--------------------------------------------------------

* **Purpose**  : The minimal interface for an object looking like a local gf function.

* **Refines**   : PureFunctionOnMesh, LazyCallable.

* **Definition** : 

* **Associated trait** : LocalGf 

LocalGfTail : the immutable local Green function's tail
------------------------------------------------------------

* **Purpose**  : The minimal interface for an object looking like the tail of a local gf function.

* **Refines**   : LazyCallable,BoostSerializable.

* **Definition** : 

+-------------------------------------------------+-------------------------------------------------------------------------------------+
| Elements                                        | Commment                                                                            |
+=================================================+=====================================================================================+
| int order_min() const int order_max() const     | The min/max order of the expansion                                                  |
+-------------------------------------------------+-------------------------------------------------------------------------------------+
| size_t size() const                             | Size ( max (0, order_max - order_min+1))                                            |
+-------------------------------------------------+-------------------------------------------------------------------------------------+
| shape_type                                      | result of shape (mini_vector<size_t,2>)                                             |
+-------------------------------------------------+-------------------------------------------------------------------------------------+
| shape_type shape() const                        | Shape of the tail                                                                   |
+-------------------------------------------------+-------------------------------------------------------------------------------------+
| mv_type operator()(size_t n) const_mv_type      | Access of the n-th order of the expansion. The non-const version throws if out of   |
| operator()(size_t n) const                      | range (order_min, order_max). The non-const version throws if out of range          |
|                                                 | (order_min, order_max)                                                              |
+-------------------------------------------------+-------------------------------------------------------------------------------------+

* **Associated trait** : LocalTail 

