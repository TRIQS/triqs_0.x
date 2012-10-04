
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

#ifndef GF_BLOC_IMLEGENDRE_qwqwegwe
#define GF_BLOC_IMLEGENDRE_qwqwegwe

#include "GF_Bloc_Base.hpp"
#include "legendre_matsubara.hpp"

typedef double Legendre_Coeffs_Type;
typedef GF_Bloc_Base<Legendre_Coeffs_Type> GF_Bloc_Legendre;

class GF_Bloc_ImLegendre : public GF_Bloc_Legendre {

public:

  GF_Bloc_ImLegendre (boost::python::object IndicesL_,
		      boost::python::object IndicesR_,
		      PyObject * Data,
		      boost::shared_ptr<MeshGF> Mesh,
		      boost::shared_ptr<TailGF> Tail);

  GF_Bloc_ImLegendre (GF_Bloc_ImLegendre const & Gin);

  const int numberLegendreCoeffs;

  PyArray<double,2> density() const;

  void determine_tail();

  void enforce_discontinuity(Array<double,2> disc);

  void enforce_discontinuity_py(python::object ob) {
    enforce_discontinuity(PyArray<double,2>(ob));
  }

 void setFromImTime (GF_Bloc_ImTime const & Gt) { legendre_matsubara_inverse(Gt,*this); }
 void setFromImFreq (GF_Bloc_ImFreq const & Gw) { legendre_matsubara_inverse(Gw,*this); }

};

#endif
