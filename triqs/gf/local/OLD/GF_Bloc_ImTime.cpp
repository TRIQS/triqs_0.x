
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

#include "GF_Bloc_ImTime.hpp"

GF_Bloc_ImTime::GF_Bloc_ImTime(boost::python::object IndicesL_,
  boost::python::object IndicesR_,
  PyObject * Data,
  boost::shared_ptr<MeshGF> Mesh,
  boost::shared_ptr<TailGF> Tail):
 GF_Bloc_Base<double>(IndicesL_,IndicesR_, Data,Mesh,Tail),
 numberTimeSlices(data.extent(2))
{
 assert(data.isStorageContiguous());
}

//------------------------------------------------------------------------------------------------

GF_Bloc_ImTime::GF_Bloc_ImTime(const GF_Bloc_ImTime & Gin):
 GF_Bloc_Base<double>(Gin),
 numberTimeSlices(Gin.numberTimeSlices)
{
 assert(data.isStorageContiguous());
}

//---------------------------------------------------------------------------

PyArray<GF_Bloc_Base<double>::element_value_type,2> GF_Bloc_ImTime::integral_tau()
{
 PyArray<GF_Bloc_Base<double>::element_value_type,2> res(N1,N2,FortranOrder);
 res=0;
 for (int i =mesh.index_min; i< mesh.index_max; i++) 
  res += data (ALL,ALL,i);
 res = res * (Beta/ mesh.len());
 return(res);
}

