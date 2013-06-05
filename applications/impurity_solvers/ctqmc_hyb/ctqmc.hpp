
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

#include <Python.h>
#include <triqs/gf/block.hpp>
#include <triqs/gf/imtime.hpp>
#include <triqs/gf/legendre.hpp>
#include <triqs/gf/imfreq.hpp>
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/mc_tools/histograms.hpp>
#include <map>
#include "configuration.hpp"

namespace triqs { namespace app { namespace impurity_solvers {

/**
   
 */
class ctqmc_hyb {

  typedef std::complex<double> SignType;

protected:

  triqs::python_tools::improved_python_dict params;
  triqs::mc_tools::HistogramBinnedMap Histograms; 
  triqs::gf::gf_view<triqs::gf::block_index,triqs::gf::gf<triqs::gf::imtime>> G_tau, F_tau, Delta_tau, OpCorrToAverage;
  triqs::gf::gf_view<triqs::gf::block_index,triqs::gf::gf<triqs::gf::legendre>> G_legendre;
  Configuration Config;
  const bool TimeAccumulation;
  const bool LegendreAccumulation;
  triqs::mc_tools::mc_generic<SignType> QMC;

public : 

  ctqmc_hyb(boost::python::object, Hloc *,
            triqs::gf::gf_view<triqs::gf::block_index,triqs::gf::gf<triqs::gf::imtime>>,
            triqs::gf::gf_view<triqs::gf::block_index,triqs::gf::gf<triqs::gf::imtime>>,
            triqs::gf::gf_view<triqs::gf::block_index,triqs::gf::gf<triqs::gf::imtime>>,
            triqs::gf::gf_view<triqs::gf::block_index,triqs::gf::gf<triqs::gf::imtime>>,
            triqs::gf::gf_view<triqs::gf::block_index,triqs::gf::gf<triqs::gf::legendre>>);
  void solve();

};

}}}

#endif
