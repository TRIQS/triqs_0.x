
/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011 by L. Boehnke, M. Ferrero, O. Parcollet
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

#include "GF_Bloc_ImLegendre.hpp"
#include <triqs/utility/legendre.hpp>

using namespace triqs::utility;

GF_Bloc_ImLegendre::GF_Bloc_ImLegendre (boost::python::object IndicesL_,
					boost::python::object IndicesR_,
					PyObject * Data,
					boost::shared_ptr<MeshGF> Mesh,
					boost::shared_ptr<TailGF> Tail):
  GF_Bloc_Legendre(IndicesL_, IndicesR_, Data, Mesh, Tail),
  numberLegendreCoeffs(data.extent(2))
{}

GF_Bloc_ImLegendre::GF_Bloc_ImLegendre (GF_Bloc_ImLegendre const & Gin):
  GF_Bloc_Legendre(Gin),
  numberLegendreCoeffs(Gin.numberLegendreCoeffs)
{}

// Compute the density from G(0^-)
PyArray<double,2> GF_Bloc_ImLegendre::density() const {
  PyArray<double,2> res(N1,N2,FortranOrder);
  for (int n1=1; n1<=N1; n1++)
    for (int n2=1; n2<=N2; n2++)
      {
	res(n1,n2)=0;
	for (int n=0; n<mesh.index_max; n++)
	  res(n1,n2)-=sqrt(2.0*n+1)*data(n1,n2,n);
      }
  res/=Beta;
  return res;
}

// Find the coefficients of the tail knowing the G_l
void GF_Bloc_ImLegendre::determine_tail() {

  tail.zero();
  for (int n1=1; n1<=N1; n1++)
    for (int n2=1; n2<=N2; n2++)
      for (int p=1; p<=tail.OrderMaxMAX; p++)
        for (int l = mesh.index_min; l<=mesh.index_max; l++)
          tail[p](n1,n2) += (legendre_t(l,p)/pow(Beta,p)) * data_const(n1,n2,l);
}

// Impose a discontinuity G(\tau=0)-G(\tau=\beta)
void GF_Bloc_ImLegendre::enforce_discontinuity(Array<double,2> disc) {
  disc.reindexSelf(TinyVector<int,2>(1,1));
  
  Array<double,1> A (Range(0,data.ubound(thirdDim)));
  for (int i=0; i<=data.ubound(thirdDim); i++)
    A(i) = legendre_t(i,1) / Beta;
  double c=sum(A*A);
  Array<double,2> step=disc.copy();
  step=disc(tensor::i,tensor::j)-sum(data(tensor::i,tensor::j,tensor::k)*A(tensor::k),tensor::k);
  data+=step(tensor::i,tensor::j)*A(tensor::k)/c;

  determine_tail();
}

