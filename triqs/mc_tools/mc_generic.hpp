
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

#ifndef TRIQS_TOOLS_MC_GENERIC_H
#define TRIQS_TOOLS_MC_GENERIC_H

#include <boost/function.hpp>
#include <math.h>
#include <triqs/utility/timer.hpp>
#include <triqs/utility/report_stream.hpp>

#include "mc_measure_set.hpp"
#include "mc_move_set.hpp"
#include "mc_basic_step.hpp"
#include "random_generator.hpp"

namespace triqs { namespace mc_tools { 

 /**
  * 
  */
 template<typename MCSignType, typename MCStepType = Step::Metropolis<MCSignType> >
 class mc_generic {

   // Couple of traits to check arguments of add_move and add_measure
   template <typename T> struct is_fine : boost::false_type {};
   template <typename T> struct is_fine<T * &&> : boost::true_type {};
   template <typename T> struct is_fine<boost::shared_ptr<T> > : boost::true_type {};
   template <typename T> struct is_fine<std::shared_ptr<T> > : boost::true_type {};

  public:

   /**
     * Constructor from a set of parameters
     */
   mc_generic(int N_Cycles, int Length_Cycle, int N_Warmup_Cycles, std::string Random_Name, int Random_Seed, int Verbosity,
              boost::function<bool()> AfterCycleDuty = boost::function<bool()>() ) :
     NCycles(N_Cycles),
     Length_MC_Cycle(Length_Cycle),
     NWarmIterations(N_Warmup_Cycles),
     RandomGenerator(Random_Name, Random_Seed),
     AllMoves(RandomGenerator),
     AllMeasures(),
     report(&std::cout, Verbosity),
     after_cycle_duty(AfterCycleDuty) {}

   /** 
     * Constructor from a dictionnary
     * \param[in] P  dictionnary parameters
     * \param[in] AfterCycleDuty  a function bool() to be called after each QMC cycle
     */
   template<typename ParameterDictType>
   mc_generic(ParameterDictType const & P, boost::function<bool()> AfterCycleDuty = boost::function<bool()>() ) :
    RandomGenerator(P.value_or_default("Random_Generator_Name",""), P.value_or_default("Random_Seed",1)),
    report(&std::cout,int(P["Verbosity"])),
    AllMoves(RandomGenerator), 
    AllMeasures(),
    Length_MC_Cycle(P["Length_Cycle"]),
    NWarmIterations(P["N_Warmup_Cycles"]),
    NCycles(P["N_Cycles"]),
    after_cycle_duty(AfterCycleDuty)
  {}

   /** 
    * Register move M with its probability of being proposed.
    * NB : the PropositionProbability needs to be >0 but does not need to be 
    * normalized. Normalization is automatically done with all the added moves
    * before starting the run
    *
      \rst
       .. warning::  
          
          The pointer is deleted automatically by the MC class at destruction. 
          Correct usage is therefore:
        
       .. code-block:: c  
      
             myMC.add_move( new myMOVE(...), Proba,Name);
    
       \endrst
    */
   template <typename MoveType>
   void add_move (MoveType * && M, std::string name, double PropositionProbability = 1.0) {
     AllMoves.add(std::move(M), name, PropositionProbability);
   }

   template <typename MoveType>
   void add_move (boost::shared_ptr<MoveType> sptr, std::string name, double PropositionProbability = 1.0) {
     AllMoves.add(sptr, name, PropositionProbability);
   }

   template <typename T>
   void add_move (T M, std::string name, double PropositionProbability = 1.0) {
     BOOST_STATIC_ASSERT_MSG(is_fine<T>::value, "Add moves with the folloing syntax: add_move(new mymove(...), name, probability) or provide a shared_ptr to a move as the first argument");
   }

   /**
    * Register the Measure M 
    *
      \rst
       .. warning::
   
          The pointer is deleted automatically by the MC class at destruction. 
          Correct usage is therefore:

       .. code-block:: c 
        
            myMC.add_measure( new myMEASURE(...), Name); 
 
      \endrst
    */
   template<typename MeasureType>
    void add_measure (MeasureType * && M, std::string name) {
      AllMeasures.insert(std::move(M), name);
    }

   template<typename MeasureType>
    void add_measure (boost::shared_ptr<MeasureType> sptr, std::string name) {
      AllMeasures.insert(sptr, name);
    }

   template<typename T>
    void add_measure (T M, std::string name) {
     BOOST_STATIC_ASSERT_MSG(is_fine<T>::value, "Add measures with the folloing syntax: add_measure(new mymeasure(...), name) or provide a shared_ptr to a measure as the first argument");
    }

   // An access to the random number generator
   random_generator RandomGenerator;

  protected:
   /**
     Reimplement to have another thermalization criterion. 
     Default is # cycles > # Warming Iterations
     */
   virtual bool thermalized() const { return (NC>= NWarmIterations);}

   /**
     Reimplement to add a convergence criterion.
     Default is false.
     It is called before each cycle and if true, the computation will stop
     */
   virtual bool converged() const { return false;}

  public:

   ///
   bool run(boost::function<bool ()> const & stop_callback) {
    assert(stop_callback);
    Timer.start();
    signe=1; done_percent = 0; 
    bool stop_it=false, finished = false;
    uint64_t NCycles_tot = NCycles+ NWarmIterations;
    report << std::endl << std::flush;
    for (NC =0; !stop_it; ++NC) {
     for (uint64_t k=1; (k<=Length_MC_Cycle); k++) { MCStepType::do_it(AllMoves,RandomGenerator,signe); }
     if (after_cycle_duty) {after_cycle_duty();}
     if (thermalized()) AllMeasures.accumulate(signe);
     // recompute fraction done
     uint64_t dp = uint64_t(floor( ( NC*100.0) / NCycles_tot));  
     if (dp>done_percent)  { done_percent=dp; report << done_percent; report<<"%; "; report <<std::flush; }
     finished = ( (NC >= NCycles_tot -1) || converged () );
     stop_it = (stop_callback() || finished);
    }
    report << std::endl << std::endl << std::flush;
    Timer.stop();
    return finished;
   }

   void collect_results(boost::mpi::communicator const & c) {
     report(2) << "Acceptance rate for all moves:" << std::endl << std::endl;
     report(2) << AllMoves.get_statistics(c) << std::endl << std::flush;
     AllMeasures.collect_results(c);
   }

  protected:

   triqs::utility::report_stream report;
   move_set<MCSignType> AllMoves;
   measure_set<MCSignType> AllMeasures;
   uint64_t Length_MC_Cycle;/// Length of one Monte-Carlo cycle between 2 measures
   uint64_t NWarmIterations,NCycles;
   triqs::utility::timer Timer;

  private: 
   boost::function<bool()> after_cycle_duty;
   MCSignType signe;
   uint64_t NC,done_percent;// NC = number of the cycle
 };


}}// end namespace
#endif

