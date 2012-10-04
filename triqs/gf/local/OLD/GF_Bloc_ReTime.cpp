
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

#include "GF_Bloc_ReTime.hpp"

// The mesh points are in the middle of the slices !
// hence tmin, tmax are not the first/last point of mesh.

 double _tmin(const MeshGF & mesh) { 
  const int N(mesh.len());
  const double Delta = real(mesh[mesh.index_max] - mesh[mesh.index_min])/(N-1);
  return real(mesh[mesh.index_min]) - Delta*0.5;
}

 double _tmax(const MeshGF & mesh) { 
  const int N(mesh.len());
  const double Delta = (real(mesh[mesh.index_max] - mesh[mesh.index_min]))/(N-1);
  return real(mesh[mesh.index_max]) + Delta*0.5;
}

GF_Bloc_ReTime::GF_Bloc_ReTime (boost::python::object IndicesL_,
				boost::python::object IndicesR_,
				PyObject * Data,
				boost::shared_ptr<MeshGF> Mesh,
				boost::shared_ptr<TailGF> Tail):
  GF_Bloc_Base<COMPLEX>(IndicesL_,IndicesR_, Data,Mesh,Tail),
  TimeMin(_tmin(mesh)), TimeMax(_tmax(mesh)),
  for_eval(data.extent(2)/(TimeMax - TimeMin))
{}

//------------------------------------------------------------------------------------------------

GF_Bloc_ReTime::GF_Bloc_ReTime(const GF_Bloc_ReTime & Gin):
  GF_Bloc_Base<COMPLEX>(Gin),
  TimeMin(Gin.TimeMin), TimeMax(Gin.TimeMax),for_eval(data.extent(2)/(TimeMax - TimeMin))
{}

//-------------------------------------------------------------------------
 COMPLEX F(COMPLEX a,double b,double Beta) {return -a/(1+exp(-Beta*b));}

PyArray<COMPLEX,2> GF_Bloc_ReTime::density() const
{
  if (!tail.is_decreasing_at_infinity()) 
    TRIQS_RUNTIME_ERROR<<" GF_Bloc_ReTime::density : Green Function is not as 1/omega or less !!!";
  
  Array<COMPLEX,2> dens_part(N1,N2,fortranArray);
  PyArray<COMPLEX,2> dens(N1,N2,FortranOrder); 
  dens_part=0;dens=0;
  
  // the < is correct : formula for  dens_part uses om+1 !
  for (int om = mesh.index_min; om<mesh.index_max; om++)
    for (int n1=1; n1<=N1;n1++) 
      for (int n2=1; n2<=N2;n2++)  
	dens_part += (data(n1,n2,om)*mesh.Bose_Fermi( real(mesh[om])) + 
		      data(n1,n2,om+1)*mesh.Bose_Fermi(real(mesh[om+1])))*real(mesh[om+1] - mesh[om])/2;
  dens_part = -imag(dens_part)/Pi;

  (Array<COMPLEX,2>)dens = dens_part; 
  return(dens);
}
 
