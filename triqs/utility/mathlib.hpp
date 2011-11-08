
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

#ifndef LIBMATH_H
#define LIBMATH_H

typedef unsigned int uint;

typedef std::complex<double> COMPLEX ;
const std::complex<double> I(0,1);
const double pi = acos(-1.0);
const double Pi = pi;
const double ZERO = 1.e-12;

#define F(op)   inline COMPLEX operator op (const COMPLEX & z, int i){return(z op double(i));} \
                inline COMPLEX operator op (int i, const COMPLEX & z){return(double(i) op z);}
F(+)
F(-)
F(/)
F(*)
#undef F


/// \f[ x^k \f]
template <class T>
T power(T x, int k)
{
  if (k==0) return(T(1));
  assert(k>0);
  T s = x;
  for (int i=1; i<k; i++) s = s * x;
  return s;
}

#endif
