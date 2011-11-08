
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

#ifndef MC_H_hewjfhwke
#define MC_H_hewjfhwke

#include "Python.h"
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/mc_tools/histograms.hpp>
#include "Configuration.hpp"
#include <map>

typedef std::complex<double> SignType;

/**
   
 */
class MC_Hybridization_Matsubara : public triqs::mc_tools::mc_generic<SignType> { 

protected:

 typedef triqs::mc_tools::mc_generic<SignType>  BaseType;
  Configuration Config;
  triqs::mc_tools::HistogramBinnedMap Histograms; 
  GF_C<GF_Bloc_ImTime> G_tau, F_tau;
  GF_C<GF_Bloc_ImLegendre> G_legendre;
  GF_C<GF_Bloc_ImFreq> Gc_w;
  const bool TimeAccumulation;
  const bool LegendreAccumulation;
  const int N_Frequencies_Accu,Freq_Fit_Start;

public : 

  MC_Hybridization_Matsubara(parameters_type const & params, size_t rank); 
  void finalize(boost::mpi::communicator const & c );

};
#endif

