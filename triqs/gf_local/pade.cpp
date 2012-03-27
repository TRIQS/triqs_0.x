
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

#include "GF_Bloc_ImFreq.hpp"
#include "GF_Bloc_ReFreq.hpp"

// This implementation is based on a Fortran code written by
// A. Poteryaev <Alexander.Poteryaev _at_ cpht.polytechnique.fr>
//
// The original algorithm is described in
// H. J. Vidberg, J. W. Serene. J. Low Temp. Phys. 29, 3-4, 179 (1977)

void define_pade_coefficients(const Array<COMPLEX,1> &z_in, const Array<COMPLEX,1> &u_in,
                              Array<COMPLEX,1> &a)
{
    int N = z_in.size();

    Array<COMPLEX,2> g(N,N);
    g = 0;
    g(0,Range::all()) = u_in;

    for(int p=1; p<N; ++p)
        for(int j=p; j<N; ++j){
            g(p,j) = (g(p-1,p-1) - g(p-1,j))/((z_in(j)-z_in(p-1))*g(p-1,j));
        }

    firstIndex i;
    a = g(i,i);
}

void define_polynomials(const Array<COMPLEX,1> &z_in, const Array<COMPLEX,1> &a,
                        Array<COMPLEX,1> &A, Array<COMPLEX,1> &B)
{
}

COMPLEX pade_substitute(const Array<COMPLEX,1> &A, const Array<COMPLEX,1> &B, COMPLEX e)
{
    // TODO
    return 0;
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

    int N1 = Gw.N1, N2 = Gw.N2;
    for (int n1=1; n1<=N1;n1++)
        for (int n2=1; n2<=N2;n2++){
            int N = Gw.mesh.len();

            Array<COMPLEX,1> u_in(N_Matsubara_Frequencies);    // Values at the Matsubara frequencies
            Array<COMPLEX,1> a(N_Matsubara_Frequencies);    // Pade coefficients

            for(int i=0; i < N_Matsubara_Frequencies; ++i){
                u_in(i) = (i < N ? Gw.data_const(n1,n2,i) : Gw.tail.eval(z_in(i))(n1,n2));
            }

            define_pade_coefficients(z_in,u_in,a);

            //Array<COMPLEX,1> A((N_Matsubara_Frequencies-1)/2 + 1)
            //Array<COMPLEX,1> B(N_Matsubara_Frequencies/2 + 1);
            //define_polynomials(z_in,a,A,B);

            int Ne = Ge.mesh.len();
            for (int i=0; i < Ne; ++i) {
                COMPLEX e = Ge.mesh[i] - I*Freq_Offset;
                //Ge.data(n1,n2,i) = pade_substitute(A,B,e);
            }
        }

    // Just copy the tail. It doesn't need to conform to the Pade approximant.
    Gw.tail = Ge.tail;
}
