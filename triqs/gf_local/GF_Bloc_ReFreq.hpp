
/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012 by M. Ferrero, O. Parcollet, I. Krivenko
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

#ifndef TRIQS_GF_BLOC_REFREQ_H
#define TRIQS_GF_BLOC_REFREQ_H

#include "GF_Bloc_Base.hpp"
#include "fourier.hpp"
#include "pade.hpp"

class GF_Bloc_ReFreq : public GF_Bloc_Base<COMPLEX> 
{
public: 
  
  GF_Bloc_ReFreq (boost::python::object IndicesL_,
		  boost::python::object IndicesR_,
		  PyObject * Data,
		  boost::shared_ptr<MeshGF> Mesh,
		  boost::shared_ptr<TailGF> Tail);
  
  GF_Bloc_ReFreq (const GF_Bloc_ReFreq & Gin);

  void setFromFourierOf(const GF_Bloc_ReTime & Gt) { fourier_direct(Gt,*this);}

  // TODO implement a more flexible way to select input Matsubara points
  void setFromPadeOf(const GF_Bloc_ImFreq & Gw, int N_Matsubara_Frequencies, double Freq_Offset);

  PyArray<COMPLEX,2> density() const;

};

#endif

