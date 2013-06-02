
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

// include first because of a namespace clash.. to be fixed...
#include "ctqmc.hpp"
#include "./iterator_python_sequence.hpp"
#include <triqs/utility/callbacks.hpp>
 
// The moves to insert and remove C, Cdagger operators
#include "move_insert.hpp"
#include "move_insert_segment.hpp"

// The move to move operators
#include "move_shift.hpp"

// The Global moves
#include "move_global.hpp"

// The measures.
#include "measure_g.hpp"
#include "measure_f.hpp"
#include "measure_average.hpp"
#include "measure_legendre.hpp"
#include "measure_correlator.hpp"

using triqs::gf::gf;
using triqs::gf::gf_view;
using triqs::gf::block_index;
using triqs::gf::imtime;
using triqs::gf::legendre;

using namespace std;
using namespace boost;


template<typename T1> 
std::string make_string( T1 const & x1) { 
 std::stringstream fs; fs <<x1; return fs.str();
}
template<typename T1,typename T2> 
std::string make_string( T1 const & x1,T2 const & x2) { 
 std::stringstream fs; fs <<x1<<x2; return fs.str();
}


namespace triqs { namespace app { namespace impurity_solvers {

//-----------------------------------------------------

ctqmc_hyb::ctqmc_hyb(boost::python::object p, Hloc * hloc,
                     gf_view<block_index,triqs::gf::gf<imtime>> gt, gf_view<block_index,triqs::gf::gf<imtime>> ft,
                     gf_view<block_index,triqs::gf::gf<imtime>> dt, gf_view<block_index,triqs::gf::gf<imtime>> opcorr,
                     gf_view<block_index,triqs::gf::gf<legendre>> gl) :
 params(p),
 G_tau(gt), F_tau(ft), Delta_tau(dt), OpCorrToAverage(opcorr), G_legendre(gl),
 Config(params, hloc, Delta_tau),
 TimeAccumulation (params["time_accumulation"]),
 LegendreAccumulation (params["legendre_accumulation"]),
 QMC(params["n_cycles"], params["length_cycle"], params["n_warmup_cycles"],
     params["random_name"], params["random_seed"], params["verbosity"])
{

 const bool UseSegmentPicture (params["use_segment_picture"]);

 // register the moves
 double p_ir = params["prob_insert_remove"];
 double p_mv = params["prob_move"];

 typedef triqs::mc_tools::move_set<SignType> move_set_type;
 move_set_type AllInserts(QMC.rng());
 move_set_type AllRemoves(QMC.rng());
 for (int a =0; a<Config.Na;++a) { 

  if (UseSegmentPicture) {
   AllInserts.add( Insert_Cdag_C_Delta_SegmentPicture ( a, Config, Histograms, QMC.rng()), make_string("Insert",a), 1.0);
   AllRemoves.add( Remove_Cdag_C_Delta_SegmentPicture ( a, Config, QMC.rng()), make_string("Remove",a), 1.0);
  }  
  else {  
   AllInserts.add(Insert_Cdag_C_Delta ( a, Config, Histograms, QMC.rng()), make_string("Insert",a), 1.0);
   AllRemoves.add(Remove_Cdag_C_Delta ( a, Config, QMC.rng()), make_string("Remove",a), 1.0);
  }
 }

 QMC.add_move(std::move(AllInserts), "INSERT", p_ir);
 QMC.add_move(std::move(AllRemoves), "REMOVE", p_ir);
 QMC.add_move(Move_C_Delta(Config, QMC.rng()), "Move C Delta", p_mv);

 // Register the Global moves
 python::list GM_List = python::extract<python::list>(params.dict()["global_moves_mapping_list"]);
 for (triqs::python_tools::IteratorOnPythonListOf3Tuples<double,python::dict,string> g(GM_List); !g.atEnd(); ++g) {
  assert (python::len(g->x2)== Config.H.N_Operators());
  // transform a python dict : name_of_operator -> name_of_operator into a  
  vector<const Hloc::Operator*> mapping(Config.H.N_Operators(), (const Hloc::Operator*)NULL);
  for (triqs::python_tools::IteratorOnPythonDict<string,string> p(g->x2); !p.atEnd(); ++p) {
   mapping[Config.H[p->key].Number] =& Config.H[p->val];
   //cout<< "MAP" << Config.H[p->key].Number<< "  "<<mapping[Config.H[p->key].Number]->Number<<endl<<
   //      Config.H[p->key].name<< "  "<<mapping[Config.H[p->key].Number]->name<<endl;
  }
  QMC.add_move(Global_Move(g->x3 , Config, QMC.rng(), mapping), "Global move", g->x1);
 }

 /*************

   Register the measures 

  ****************/

 for (int a =0; a<Config.Na;++a) { 
   if (LegendreAccumulation) {
     QMC.add_measure(Measure_G_Legendre(Config, a, G_legendre[a]), make_string("G Legendre ",a));
   } else if (TimeAccumulation) {
     QMC.add_measure(Measure_G_tau(Config, a, G_tau[a] ), make_string("G(tau) ",a));
   } else {
     assert(0);
   }
 }


 // register the measure of F
 if (bool(params["use_f"]))
   for (int a =0; a<Config.Na;++a) 
     QMC.add_measure(Measure_F_tau(Config, a, F_tau[a] ), make_string("F(tau) ",a));

 // register the measures of the average of some operators
 python::dict opAv_results = python::extract<python::dict>(params.dict()["measured_operators_results"]);
 python::list opAv_List = python::extract<python::list>(params.dict()["operators_av_list"]);
 for (triqs::python_tools::IteratorOnPythonList<string> g(opAv_List); !g.atEnd(); ++g) {
  QMC.add_measure(Measure_OpAv(*g, Config, opAv_results), *g);
 }

 // register the measures for the time correlators:
 python::list opCorr_List = python::extract<python::list>(params.dict()["opcorr_av_list"]);
 int a = 0;
 for (triqs::python_tools::IteratorOnPythonList<string> g(opCorr_List); !g.atEnd(); ++g, ++a) {
  string str1(*g);
  str1+= "OpCorr";
  QMC.add_measure(Measure_OpCorr(str1, *g, Config, OpCorrToAverage[a], OpCorrToAverage[a].mesh().size()), str1);
 }

}


//********************************************************

void ctqmc_hyb::solve() {

  // communicate on the world = all nodes
  boost::mpi::communicator c;

  // run!! The empty configuration has sign = 1
  QMC.start(1.0, triqs::utility::clock_callback(params["max_time"]));
  QMC.collect_results(c);

}

}}}
