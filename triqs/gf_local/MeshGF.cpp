
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

#include "MeshGF.hpp"
using namespace boost;
using namespace std;
using python::extract;
using python::object;
using python::dict;

MeshGF::MeshGF(const Type_GF & typeGF_, const Statistic_GF & stat, double Beta_, PyObject * MeshArray):
  arr(MeshArray),//view
  Beta(Beta_),
  Statistic(stat),
  typeGF(typeGF_),
  index_min(0),index_max(arr.extent(0)-1)
{ 
  iter = index_min;
}

MeshGF::MeshGF(const MeshGF & X) :
  arr(X.arr.as_PyObjectPtrBorrowedRef()), //no copy, array is never changed in the class
  Beta(X.Beta),
  Statistic(X.Statistic),
  typeGF(X.typeGF),
  index_min(X.index_min),
  index_max(X.index_max)
{
  iter = index_min;
}

MeshGF::MeshGF(const MeshGF & X, const Statistic_GF & stat) :
  arr(X.arr.as_PyObjectPtrBorrowedRef()), //no copy, array is never changed in the class
  Beta(X.Beta),
  Statistic(stat),
  typeGF(X.typeGF),
  index_min(X.index_min),
  index_max(X.index_max)
{
  iter = index_min;
}

double MeshGF::Bose_Fermi(double x) const {
    const double coupure=1e-12,cutoff=100;
    double  y=Beta*x;
    if (Statistic==Boson) {
      if ((abs(y)<coupure)||(y>cutoff)) return(0);
      if ((-y)>cutoff) return(-1);
      return( 1/(exp(y) - 1)); }
    else
      return (1/(exp(y)+1));
  }

object MeshGF::next(){
  if (iter>index_max) {
    PyErr_SetString(PyExc_StopIteration, "");
    python::throw_error_already_set();
  }
  object r( (typeGF == Imaginary_Frequency ? object(I*arr(iter)) : object(arr(iter))) ); 
  iter++; return r;
}

object MeshGF::__reduce__() const { 
  return python::make_tuple( python::import("pytriqs.base.utility.myUtils").attr("call_factory_from_dict"),
			     python::make_tuple(object(*this).attr("__class__"),__reduce_to_dict__()));

}

object MeshGF::__reduce_to_dict__() const {
  dict d; // construct a new dictionnary
  d["array"] = object(arr);
  d["TypeGF"] = object(typeGF).attr("name");
  d["Statistic"] = object(Statistic).attr("name");
  d["Beta"]= object(Beta);
  return d;
}

object MeshGF::__factory_from_dict__(const object & dic){
  dict d(dic);
  object arr(d["array"]);
  Type_GF t1(Imaginary_Frequency); object cl1= object(t1).attr("__class__");
  string s1 = extract<string>(d["TypeGF"]);
  Type_GF t( extract<Type_GF>( cl1.attr(s1.c_str()) ));
  Statistic_GF t2(Fermion); object cl2= object(t2).attr("__class__");
  string s2 = extract<string>(d["Statistic"]);
  Statistic_GF s( extract<Statistic_GF>( cl2.attr(s2.c_str()) ));
  double beta = extract<double>(d["Beta"]);
  return object(MeshGF(t,s,beta,arr.ptr()));
}

bool MeshGF::operator==(const MeshGF & othermesh) const {
  return  (abs(Beta - othermesh.Beta)<ZERO) &&
    (Statistic ==othermesh.Statistic) &&
    (typeGF ==othermesh.typeGF) &&
    ((index_min == othermesh.index_min) && (index_max == othermesh.index_max)) &&
    (max(abs(arr- othermesh.arr))<ZERO);
}  

