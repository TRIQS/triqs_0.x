
/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012 by M. Ferrero, O. Parcollet
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

#ifndef TRIQS_GF_BLOC_PADE
#define TRIQS_GF_BLOC_PADE

#include "GF_Bloc_ImFreq.hpp"
#include "GF_Bloc_ReFreq.hpp"

namespace gmp{
#include <gmpxx.h>
}

typedef std::complex<gmp::mpf_class> MP_COMPLEX;

class Pade_approximant {

    Array<COMPLEX,1> z_in; // Input complex frequency points 
    Array<COMPLEX,1> a; // Pade coefficients
    
    public:
    
    Pade_approximant(const Array<COMPLEX,1> &z_in, const Array<COMPLEX,1> &u_in);
    COMPLEX operator()(COMPLEX e) const; // Calculate the value of the approximant at a given point
    
    static const int GMP_default_prec = 256;    // Precision of GMP floats to use during a Pade coefficients calculation.
};

void pade (GF_Bloc_ImFreq const & Gw, GF_Bloc_ReFreq & Ge, int N_Matsubara_Frequencies, double Freq_Offset);

#endif
