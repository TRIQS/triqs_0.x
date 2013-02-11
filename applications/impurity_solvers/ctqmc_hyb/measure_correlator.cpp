
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

#include "measure_z.hpp"
#include "measure_correlator.hpp"

void Measure_OpCorr::accumulate(std::complex<double> signe) { 

   const double s(real(signe)); // not elegant !
   BaseType::accumulate(s);
  double tau0 = 0.00, taueval, tauold=0.0, taunew; 
  int indlo=0, indhi, ii;
  
  for (Configuration::OP_REF op = Config.DT.OpRef_begin(); ! op.atEnd(); ++op) {
    // set time of evaluation in the middle of two operators:
    taunew = op->tau;
    taueval = (taunew + tauold) / 2.0;
    // now calculate the traces:
    Config.DT.insertTwoOperators(taueval,Config.H[opName],tau0,Config.H[opName]);
    double val = Config.DT.ratioNewTrace_OldTrace();
    Config.DT.undo_insertTwoOperators();

    // bin the result:
    indhi = int(floor(abs(taunew/deltatau)));
    if (indlo!=indhi) {
      // put data into Op_res:
      Op_res_bin(0, (indlo+0.5)*deltatau, 0, 0.0, ( (indlo+1)*deltatau - tauold) * val * s );
      Op_res_bin(0, (indhi+0.5)*deltatau, 0, 0.0, ( taunew - indhi*deltatau) * val * s );
      for (ii=indlo+1; ii<indhi; ii++) {
        Op_res_bin(0, (ii+0.5)*deltatau, 0, 0.0, deltatau * val * s);
      }
    }
    else {
      Op_res_bin(0, (indlo+0.5)*deltatau, 0, 0.0, (taunew - tauold) * val * s);
    }

    tauold = taunew;
    indlo = indhi;
  }
  
  // Add last segment:
  taunew = Config.Beta;
  taueval = (taunew + tauold) / 2.0;
  // now calculate the traces:
  Config.DT.insertTwoOperators(taueval,Config.H[opName],tau0,Config.H[opName]);
  double val = Config.DT.ratioNewTrace_OldTrace();
  Config.DT.undo_insertTwoOperators();
  // bin the result:
  indhi = N_timeslices-1;
  if (indlo!=indhi) {
    // put data into Op_Res:
    Op_res_bin(0, (indlo+0.5)*deltatau, 0, 0.0, ( (indlo+1)*deltatau - tauold) * val * s );
    Op_res_bin(0, (indhi+0.5)*deltatau, 0, 0.0, ( taunew - indhi*deltatau) * val * s );
    for (ii=indlo+1; ii<indhi; ii++) {
      Op_res_bin(0, (ii+0.5)*deltatau, 0, 0.0, deltatau * val * s);
    }
  }
  else {
    Op_res_bin(0, (indlo+0.5)*deltatau, 0, 0.0, (taunew - tauold) * val * s);
  }
}

