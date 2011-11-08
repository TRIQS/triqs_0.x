
/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011 by M. Ferrero, O. Parcollet
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#ifndef TRIQS_GF_BLOC_BASE
#define TRIQS_GF_BLOC_BASE

#include "MeshGF.hpp"
#include "TailGF.hpp"
#include <boost/shared_ptr.hpp>
#include <triqs/utility/exceptions.hpp>

template <typename DataType>
class GF_Bloc_Base {
 public: 
  /*
     - IndicesL,R : indices of the left/right index.
     - Data : anything that numpy can transform into a 3d complex array of dim (N1,N2, Nmax)
     - Mesh : a meshGF object viewed in python
     - Tail : a TailGF object viewed in python
     */
  GF_Bloc_Base (python::object IndicesL_,
    python::object IndicesR_,
    PyObject * Data,
    boost::shared_ptr<MeshGF> Mesh,
    boost::shared_ptr<TailGF> Tail);

  ~GF_Bloc_Base();

  // Returns a VIEW of the same function !
  GF_Bloc_Base (const GF_Bloc_Base<DataType> & Gin);

  // WARNING : indices are NOT COPIED in these routines. They are CONST
  void operator= (const GF_Bloc_Base<DataType> & Gin);

  const python::list IndicesL, IndicesR;

  boost::shared_ptr<MeshGF> mesh_ptr;
  boost::shared_ptr<TailGF> tail_ptr;
  MeshGF & mesh; // the mesh of the Green function
  TailGF & tail; // the large omega expansion

  const double Beta; // Copy of the mesh data for convenience
  const Statistic_GF Statistic; // Copy of the mesh data for convenience
  const int N1, N2;

  typedef DataType element_value_type;

 private:
  PyArray<DataType,3> * data_ptr;

 public : // not safe : replace by a view which can not be resized....
  // but still one wishes to have access to the data.
  // solved when replacing blitz++ by new array class...
  //protected : 
  PyArray<DataType,3> & data;

 public :

 // PyArray<DataType,3> _data()      { return data;}
  const PyArray<DataType,3> & data_const;

  // access to data for python interface
  python::object data_as_numpy() const { return data.as_BoostObject();}

  void zero() {data =0.0;tail.zero();}

  // reduce all the arrays on the master. Requires not distributed
  void MPI_reduce_sum_onsite();

  // Bcast the data over the nodes and set the tail on the node.
  void MPI_bcast();

  /// Save the Green function in i omega_n (as 2 columns).
  void save(string file,  bool accumulate=false) const;
  /// Load the GF
  void load(string file);

  void operator *= (DataType alpha) {data *= alpha; tail *=alpha;}
  void operator /= (DataType alpha) {data /= alpha; tail /=alpha;}

};

 template<typename T1, typename T2>
void check_have_same_structure(GF_Bloc_Base<T1> const & G1, GF_Bloc_Base<T2> const & G2, bool mesh_check=true, bool check_NN= true) 
{ 
 if (mesh_check) G1.mesh.check_is_same(G2.mesh);
 if (check_NN && (G1.N1!= G2.N1))  TRIQS_RUNTIME_ERROR<<"The first dimension of the two Green functions is not the same";
 if (check_NN && (G1.N2!= G2.N2))  TRIQS_RUNTIME_ERROR<<"The second dimension of the two Green functions is not the same";
 if (abs(G1.Beta - G2.Beta)> ZERO) TRIQS_RUNTIME_ERROR<<"The two Green functions have not the same Beta"<<"G1 :"<<G1.Beta<<" G2 : "<<G2.Beta;
}


#endif

