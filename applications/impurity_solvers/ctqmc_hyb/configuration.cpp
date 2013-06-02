
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

#include "ctqmc.hpp"
#include "./iterator_python_sequence.hpp"

using namespace std;
using namespace boost;
using namespace ctqmc_utils;
using triqs::gf::gf;
using triqs::gf::gf_view;
using triqs::gf::block_index;
using triqs::gf::imtime;

Configuration::Configuration(triqs::python_tools::improved_python_dict params, Hloc * hloc, gf_view<block_index,gf<imtime>> & dt):
  Delta_tau(dt),
  H(*hloc),
  Beta(dt[0].domain().beta),
  DT(H,Beta), 
  Na(python::len(params.dict()["c_cdag_ops"])),
  dets(Na,(DET_TYPE *) NULL),
  COps(Na),CdagOps(Na),
  Delta_tau_proxy(Na,(Delta_Proxy*) NULL), 
  info (H.N_Operators()),
  RecordStatisticConfigurations(params["record_stat"]),
  CurrentSign(1),
  OldSign(1)
{
  // c_cdag_ops is expected to be a list of list of NAMES of operators (C,C_dag)
  // prepared in the python solver
  python::list BOL = python::extract<python::list>(params.dict()["c_cdag_ops"]);
  int a=0;
  for (triqs::python_tools::IteratorOnPythonList<python::list> g(BOL); !g.atEnd(); ++g, ++a) {
    int alpha=0;
    for (triqs::python_tools::IteratorOnPythonListOf2Tuples<string,string> p(*g); !p.atEnd(); ++p, ++alpha) {
      COps[a].push_back(&H[p->x1]);
      info[COps[a][alpha]->Number].a = a;
      info[COps[a][alpha]->Number].alpha = alpha;
      info[COps[a][alpha]->Number].dagger = false;

      CdagOps[a].push_back(&H[p->x2]);
      info[CdagOps[a][alpha]->Number].a = a;
      info[CdagOps[a][alpha]->Number].alpha = alpha;
      info[CdagOps[a][alpha]->Number].dagger = true;
    }
  }
    
  // construct the determinants
  int Nmax = params["n_max_matrix"];
  for (int a =0; a<Na;++a) { 
    Delta_tau_proxy[a] = new Delta_Proxy(Delta_tau[a],info);
    dets[a] = new DET_TYPE(*Delta_tau_proxy[a],Nmax);
  }
  
}

//********************************************************

Configuration::~Configuration() { 
  for (uint i =0; i<dets.size();i++) {delete dets[i];delete Delta_tau_proxy[i];} 
}
 
//********************************************************

Configuration::O_Odag_Insertions_type & Configuration::create_O_Odag_Insertion (const Hloc::Operator & Op1, const Hloc::Operator & Op2) { 
  map<std::pair<string,string>, O_Odag_Insertions_type >::iterator ret; bool ok;
  std::tie(ret,ok) = O_Odag_Insertions.insert(std::make_pair( std::make_pair(Op1.name,Op2.name), Configuration::O_Odag_Insertions_type()));
  if (!ok) TRIQS_RUNTIME_ERROR<<"Double insert of O_Odag_insertion !!!!";
  O_Odag_Insertions_type & RES((*ret).second);
  RES.clear();
  this->info[Op1.Number].InsertionTablePtr =  &( RES.Insertions1);
  this->info[Op2.Number].InsertionTablePtr =  &( RES.Insertions2);
  return RES;
}

//********************************************************

// Compute the sign of the config from scratch
void Configuration::update_Sign() {

  int s =0;
  vector<int> n_op_with_a(Na,0), n_op_with_dag_a(Na,0);

  // In this first part we compute the sign to bring the configuration to
  // d^_1 d^_1 d^_1 ... d_1 d_1 d_1   d^_2 d^_2 ... d_2 d_2   ...   d^_n .. d_n

  // loop over the operators "op" in the trace (right to left)
  for (Configuration::OP_REF op = DT.OpRef_begin(); ! op.atEnd(); ++op) {

    // the info of the operator in "op"
    BlockInfo op_info(info[op->Op->Number]);

    // Check
    if ((op_info.a == -1) && (op->Op->Statistic != Hloc::Operator::Bosonic))
      TRIQS_RUNTIME_ERROR<<"Sign Computation : Operators must be fundamental or bosonic !!";

    //DEBUG. check fundamental ops.
    if (op_info.a != -1) assert (op->Op->Statistic == Hloc::Operator::Fermionic);

    // Skip the bosonic operator in the sign computation
    if (op->Op->Statistic == Hloc::Operator::Bosonic) continue;

    // how many operators with a smaller a than "op" are there on the right of "op"?
    for (int a = 0; a < op_info.a; a++) s += n_op_with_a[a];
    n_op_with_a[op_info.a]++;

    // if "op" is not a dagger how many operators of the same a but with a dagger are there on his right?
    if (!op_info.dagger)
      s += n_op_with_dag_a[op_info.a];
    else
      n_op_with_dag_a[op_info.a]++;
  }

  // Now we compute the sign to bring the configuration to
  // d_1 d^_1 d_1 d^_1 ... d_1 d^_1   ...   d_n d^_n ... d_n d^_n

  for (int a = 0; a < Na; a++) {
    int n = dets[a]->size();
    s +=  n*(n+1)/2;
  }

  OldSign = CurrentSign;
  CurrentSign = (s%2==0 ? 1 : -1);

}
