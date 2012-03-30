
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

#ifndef TRIQS_BASE_GF_BLOC_IMFREQ_H
#define TRIQS_BASE_GF_BLOC_IMFREQ_H

#include "GF_Bloc_Base.hpp"
#include "fourier.hpp"
#include "legendre_matsubara.hpp"

class GF_Bloc_ImFreq : public GF_Bloc_Base<COMPLEX> {

public: 
  
  GF_Bloc_ImFreq (boost::python::object IndicesL_,
		  boost::python::object IndicesR_,
		  PyObject * Data,
		  boost::shared_ptr<MeshGF> Mesh,
		  boost::shared_ptr<TailGF> Tail);
  
  GF_Bloc_ImFreq (const GF_Bloc_ImFreq & Gin);

  /// Computes the density \f$ G(\tau = 0^-) \f$
  PyArray<COMPLEX,2> density() const;

  void setFromFourierOf(const GF_Bloc_ImTime & Gt, bool time_mesh_starts_at_half_bin = true) { fourier_direct(Gt,*this,time_mesh_starts_at_half_bin); }
  void setFromLegendre(const GF_Bloc_ImLegendre & Gl) { legendre_matsubara_direct(Gl,*this); }
 
};

#endif
















