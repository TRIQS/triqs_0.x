
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

#ifndef TRIQS_GF_BLOC_FOURIER
#define TRIQS_GF_BLOC_FOURIER

class GF_Bloc_ImTime;
class GF_Bloc_ImFreq;
class GF_Bloc_ReTime;
class GF_Bloc_ReFreq;

void fourier_direct  (GF_Bloc_ImTime const & Gt, GF_Bloc_ImFreq & Gw, bool time_mesh_starts_at_half_bin);
void fourier_inverse (GF_Bloc_ImFreq const & Gw, GF_Bloc_ImTime & Gt, bool time_mesh_starts_at_half_bin);

void fourier_direct  (GF_Bloc_ReTime const & Gt, GF_Bloc_ReFreq & Gw);
void fourier_inverse (GF_Bloc_ReFreq const & Gw, GF_Bloc_ReTime & Gt );

#endif




