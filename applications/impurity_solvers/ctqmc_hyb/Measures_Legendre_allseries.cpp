
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

#include "Measures_Legendre_allseries.hpp"
#include <triqs/utility/legendre.hpp>

using namespace triqs::utility;

void Measure_G_Legendre_all::accumulate(COMPLEX signe) {
 
  const double s(real(signe)); // not elegant !
  BaseType::accumulate(s);

  legendre_generator Tn = legendre_generator();
  double delta_tau;
  short st;

  data_part()=0;
  for (Configuration::DET_TYPE::C_Cdagger_M_iterator p(*conf.dets[a_level]); !p.atEnd(); ++p) {

      st=1;
      delta_tau=p.C()->tau-p.Cdagger()->tau;
      if (delta_tau<0.0) {
        delta_tau+=conf.Beta;
        st=-1;
      }

      const double x=2*delta_tau/conf.Beta-1;
      Tn.reset(x);
      const double temp=s*st*p.M();
      for (int n=0; n<Gl.numberLegendreCoeffs; n++) {
        double a =  temp*Tn.next();
        data_part(conf.info[p.C()->Op->Number].alpha, conf.info[p.Cdagger()->Op->Number].alpha, n) += a;
      }

  }
  data +=data_part;
  data_stack << data_part;
  
}

