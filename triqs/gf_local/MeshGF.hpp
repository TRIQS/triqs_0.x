
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

#ifndef TRIQS_BASE_GF_LOCAL_MESH_H 
#define TRIQS_BASE_GF_LOCAL_MESH_H

#include "Python.h"
#include <boost/python.hpp>
#ifndef NO_USE_BLITZ
  #include "triqs/deprecated/blitzext/PyArray.hpp"
#endif
#include <triqs/utility/mathlib.hpp>
#include <triqs/utility/exceptions.hpp>

enum Type_GF {Imaginary_Frequency,Real_Frequency,Imaginary_Time,Real_Time,Imaginary_Legendre};

// Explicitly give the values consistent with 'XOR multiplication' of statistics.
enum Statistic_GF {Boson=0,Fermion=1};

/*
  Contains : 
     - the mesh array
     - beta
     - type of GF and statistic
*/
class MeshGF { 

  const PyArray<double,1> arr; // The underlying array

public :

  // MeshArray is a numpy array containing the mesh. A view of it will be taken.
  MeshGF(const Type_GF & typeGF_, const Statistic_GF & stat, double Beta_, PyObject * MeshArray);

  // Another view on the same array
  MeshGF(const MeshGF & X);
  
  const double Beta; 
  Statistic_GF Statistic;
  const Type_GF typeGF;     
  const int index_min, index_max; 

  inline int len() const {return index_max - index_min + 1;} // length
  
  // The Mesh is an iterator
  MeshGF __iter__(){ return *this;}  
  boost::python::object next();

  // access to the value // COMPLEX OR REAL : to be discussed 
  std::complex<double> operator[](int i) const { assert( (i<=index_max) && (i>=index_min)); 
    return (typeGF == Imaginary_Frequency ? I*arr(i) : arr(i));}

  bool operator==(const MeshGF &other) const;

  double Bose_Fermi(double x) const;

  inline void check_is_same (const MeshGF & othermesh) { 
    if (!(*this == othermesh)) TRIQS_RUNTIME_ERROR<<"Mesh are not the same";}

  boost::python::object __reduce__() const;
  boost::python::object __reduce_to_dict__() const;
  static boost::python::object __factory_from_dict__(const boost::python::object & dic);
  
private:
  template<class T> void operator = (const T& x); //forbidden
  int iter,stop;
};

#endif
