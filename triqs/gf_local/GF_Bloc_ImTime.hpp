
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

#ifndef TRIQS_BASE_GF_BLOC_TIME_H
#define TRIQS_BASE_GF_BLOC_TIME_H

#include "GF_Bloc_Base.hpp"
#include "fourier.hpp"
#include "legendre_matsubara.hpp"

class GF_Bloc_ImTime : public GF_Bloc_Base<double> {
 void operator= (const GF_Bloc_ImTime & Gin); // not implemented to be fixed. there are const here....
 public: 
 GF_Bloc_ImTime (boost::python::object IndicesL_,
   boost::python::object IndicesR_,
   PyObject * Data,
   boost::shared_ptr<MeshGF> Mesh,
   boost::shared_ptr<TailGF> Tail);

 GF_Bloc_ImTime (const GF_Bloc_ImTime & Gin);
 
 const int numberTimeSlices;

 void setFromInverseFourierOf(const GF_Bloc_ImFreq & Gw, bool time_mesh_starts_at_half_bin = true) { fourier_inverse(Gw, *this,time_mesh_starts_at_half_bin);}
 void setFromLegendre (GF_Bloc_ImLegendre const & Gl) { legendre_matsubara_direct(Gl,*this); }

 PyArray<GF_Bloc_Base<double>::element_value_type,2> integral_tau();

};
#endif

