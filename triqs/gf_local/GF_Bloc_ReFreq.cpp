
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

#include "GF_Bloc_ReFreq.hpp"
#include "GF_Bloc_ReTime.hpp"

GF_Bloc_ReFreq::GF_Bloc_ReFreq (boost::python::object IndicesL_,
				boost::python::object IndicesR_,
				PyObject * Data,
				boost::shared_ptr<MeshGF> Mesh,
				boost::shared_ptr<TailGF> Tail):
  GF_Bloc_Base<COMPLEX>(IndicesL_,IndicesR_, Data,Mesh,Tail)
{}
 
//------------------------------------------------------------------------------------------------

GF_Bloc_ReFreq::GF_Bloc_ReFreq(const GF_Bloc_ReFreq & Gin):
  GF_Bloc_Base<COMPLEX>(Gin)
{}

//-------------------------------------------------------------------------
inline COMPLEX F(COMPLEX a,double b,double Beta) {return -a/(1+exp(-Beta*b));}

PyArray<COMPLEX,2> GF_Bloc_ReFreq::density() const
{
  if (!tail.is_decreasing_at_infinity()) 
    TRIQS_RUNTIME_ERROR<<" GF_Bloc_ReFreq::density : Green Function is not as 1/omega or less !!!";
  
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
 
void GF_Bloc_ReFreq::setFromPadeOf(const GF_Bloc_ImFreq & Gw, int N_Matsubara_Frequencies, double Freq_Offset)
{
    if(Freq_Offset < 0) TRIQS_RUNTIME_ERROR << "Frequency offset must be non-negative (Freq_Offset = " << Freq_Offset << ")";

    pade(Gw,*this,N_Matsubara_Frequencies,Freq_Offset);
}