
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

#include "pade.hpp"

// This implementation is based on a Fortran code written by
// A. Poteryaev <Alexander.Poteryaev _at_ cpht.polytechnique.fr>
//
// The original algorithm is described in
// H. J. Vidberg, J. W. Serene. J. Low Temp. Phys. 29, 3-4, 179 (1977)

Pade_approximant::Pade_approximant(const Array<COMPLEX,1> & z_in, const Array<COMPLEX,1> & u_in) :
    z_in(z_in)
{
    int N = z_in.size();
    a.resize(N);
  
    // Change the default precision of GMP floats.
    // MF: I have put an unsigned long instead of gmp::mp_bitcnt_t here so that it's compatible with GMP >= 3.0
    unsigned long old_prec = gmp::mpf_get_default_prec();
    gmp::mpf_set_default_prec(GMP_default_prec);  // How do we determine it?
 
    Array<MP_COMPLEX,2> g(N,N);
    g = MP_COMPLEX(.0);
    g(0,Range::all()) = cast<MP_COMPLEX>(u_in);
 
    MP_COMPLEX MP_1(1.0);
    
    for(int p=1; p<N; ++p)
        for(int j=p; j<N; ++j){
            MP_COMPLEX x(g(p-1,p-1)/g(p-1,j) - MP_1);
            MP_COMPLEX y(z_in(j)-z_in(p-1));
            g(p,j) = x/y;
        }

    for(int j=0; j<N; ++j){
        MP_COMPLEX gj(g(j,j));
        a(j) = COMPLEX(real(gj).get_d(),imag(gj).get_d());
    }
    
    // Restore the precision.
    gmp::mpf_set_default_prec(old_prec);
}

COMPLEX Pade_approximant::operator()(COMPLEX e) const
{
    COMPLEX A1(0);
    COMPLEX A2 = a(0);
    COMPLEX B1(1.0), B2(1.0);

    int N = a.size();
    for(int i=0; i<=N-2; ++i){
        COMPLEX Anew = A2 + (e - z_in(i))*a(i+1)*A1;
        COMPLEX Bnew = B2 + (e - z_in(i))*a(i+1)*B1;
        A1 = A2; A2 = Anew;
        B1 = B2; B2 = Bnew;
    }

    return A2/B2;
}

void pade(GF_Bloc_ImFreq const & Gw, GF_Bloc_ReFreq & Ge, int N_Matsubara_Frequencies, double Freq_Offset)
{
    check_have_same_structure (Gw,Ge,false,true);
    assert (Gw.mesh.index_min==0);
    assert (Ge.mesh.index_min==0);

    double Beta = Gw.Beta;
    double omegaShift = (Gw.Statistic==Fermion ? 1 : 0);

    Array<COMPLEX,1> z_in(N_Matsubara_Frequencies);
    firstIndex i;
    z_in = I*Pi/Beta*(2*i+omegaShift);

    // Just copy the tail. It doesn't need to conform to the Pade approximant.
    Gw.tail = Ge.tail;

    int N1 = Gw.N1, N2 = Gw.N2;
    for (int n1=1; n1<=N1;n1++)
        for (int n2=1; n2<=N2;n2++){
            int N = Gw.mesh.len();

            Array<COMPLEX,1> u_in(N_Matsubara_Frequencies);     // Values at the Matsubara frequencies
            Array<COMPLEX,1> a(N_Matsubara_Frequencies);        // Pade coefficients

            for(int i=0; i < N_Matsubara_Frequencies; ++i){
                u_in(i) = (i < N ? Gw.data_const(n1,n2,i) : Gw.tail.eval(z_in(i))(n1,n2));
            }

            Pade_approximant PA(z_in,u_in);

            int Ne = Ge.mesh.len();
            Ge.zero();
            for (int i=0; i < Ne; ++i) {
                COMPLEX e = Ge.mesh[i] + I*Freq_Offset;
                Ge.data(n1,n2,i) = PA(e);
            }
        }
}
