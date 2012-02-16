
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

#ifndef TRIQS_BASE_GF_BLOC_RETIME_H
#define TRIQS_BASE_GF_BLOC_RETIME_H

#include "GF_Bloc_Base.hpp"
#include "fourier.hpp"

class GF_Bloc_ReTime :  public GF_Bloc_Base<COMPLEX>
{
 public: 

  GF_Bloc_ReTime (
    boost::python::object IndicesL_,
    boost::python::object IndicesR_,
    PyObject * Data,
    boost::shared_ptr<MeshGF> Mesh,
    boost::shared_ptr<TailGF> Tail);

  GF_Bloc_ReTime (const GF_Bloc_ReTime & Gin);

  // Linear evaluation 
  COMPLEX operator() (double t, int n1, int n2) const { 
   int i = int(floor((t-TimeMin)*for_eval -0.5 )); // cf .py, the mesh is centered
   if ((i<data.lbound(2)) ||(i>=data.ubound(2))) return 0;
   double t0 = TimeMin + (i+0.5)/for_eval;
   double x = (t-t0)*for_eval; assert(x>=0); assert(x<=1);
   return data(n1,n2,i) * (1-x) + x*data(n1,n2,i+1);
  }

  /// Compute the inverse Fourier transform from Gw 
  void setFromInverseFourierOf(const GF_Bloc_ReFreq & Gw) { fourier_inverse(Gw,*this);}

  /// Computes the density \f$ G(\tau = 0^-) \f$
  PyArray<COMPLEX,2> density() const;

  /// Time Window
  const double TimeMin,TimeMax;
 private:
  const double for_eval;
};

#endif
















