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
#ifndef TRIQS_GF_LOCAL_FOURIER_REAL_H 
#define TRIQS_GF_LOCAL_FOURIER_REAL_H

#include <triqs/gf/descriptors/freq.hpp> 
#include <triqs/gf/descriptors/one_time.hpp> 

namespace triqs { namespace gf { 

gf<freq> fourier_direct (gf<one_time> const & gt);
gf<one_time> fourier_inverse (gf<freq> const & gw);

}}
#endif


