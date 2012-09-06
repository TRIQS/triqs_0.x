
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

#ifndef TRIQS_TOOLS_MC_MOVE_SET2_H
#define TRIQS_TOOLS_MC_MOVE_SET2_H
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/function.hpp>
#include <boost/mpi.hpp>
#include <boost/concept_check.hpp>
#include <boost/utility/enable_if.hpp>
#include "polymorphic_random_generator.hpp"
#include <triqs/utility/report_stream.hpp>
#include <triqs/utility/exceptions.hpp>

namespace triqs { namespace mc_tools { 
 namespace mpi=boost::mpi; namespace BLL = boost::lambda;

 template <class X, typename Y> struct IsMove { // a concept check for the moves at compile time
  BOOST_CONCEPT_USAGE(IsMove)
  {
   Y r = i.Try();      
   r = i.Accept();
   i.Reject();
  }
  private: X i;
 };

 template<typename MCSignType> class move_set;

 template<typename MCSignType> 
  class move {
   boost::shared_ptr<void> move_impl;
   uint64_t NProposed,NAccepted; // Statistics
   move_set<MCSignType> * mset_ptr; // if hte move is a move_set, we keep a handle for print reports ...
   // two little helper functions for the constructor : pointer is not NULL iif we give a moveset, otherwise it is NULL.
   static move_set<MCSignType> * _mk_ptr( move_set<MCSignType> * x) { return x;}
   template<class T> static move_set<MCSignType> * _mk_ptr( T * x) { return NULL;}

   public:
   boost::function<MCSignType()> Try_, Accept_;
   boost::function<void()> Reject;

   template<typename MoveType>
    move( MoveType * move_ptr) : 
     move_impl(move_ptr),NProposed(0),NAccepted(0),
     mset_ptr( _mk_ptr(move_ptr)), 
     Try_(BLL::bind(&MoveType::Try,move_ptr)),
     Accept_(BLL::bind(&MoveType::Accept,move_ptr)),
     Reject(BLL::bind(&MoveType::Reject,move_ptr))
   {
    BOOST_CONCEPT_ASSERT((IsMove<MoveType,MCSignType>));
   }

   template<typename MoveType>
    move(boost::shared_ptr<MoveType> sptr) :
     move_impl(sptr),NProposed(0),NAccepted(0),
     mset_ptr( _mk_ptr(sptr.get())),
     Try_(BLL::bind(&MoveType::Try,sptr.get())),
     Accept_(BLL::bind(&MoveType::Accept,sptr.get())),
     Reject(BLL::bind(&MoveType::Reject,sptr.get()))
   {
    BOOST_CONCEPT_ASSERT((IsMove<MoveType,MCSignType>));
   }

   MCSignType Try(){ NProposed++; return Try_();}
   MCSignType Accept() { NAccepted++; return  Accept_(); }

   double acceptance_probability(mpi::communicator const & c) const { 
    uint64_t nacc_tot=0, nprop_tot=1;
    mpi::reduce(c, NAccepted, nacc_tot, std::plus<uint64_t>(), 0);
    mpi::reduce(c, NProposed, nprop_tot, std::plus<uint64_t>(), 0);
    return nacc_tot/static_cast<double>(nprop_tot);
   }

   void print (triqs::utility::report_stream & report, mpi::communicator const & c, std::string name, std::string decal) const {
    report<< decal <<"Acceptance probability of move : "<<name<<"  :  "<<acceptance_probability(c)<<std::endl;
    if (mset_ptr) mset_ptr->print (report, c, name, decal);
   }
  };// class move  

 //--------------------------------------------------------------------

 /// A vector of (moves, PropositionProbability), which is also a move itself
 template<typename MCSignType>
  class move_set : std::vector<move<MCSignType> > { 
   std::vector<std::string> names_;
   move<MCSignType> * current;
   size_t current_move_number;
   polymorphic_random_generator & RNG;
   std::vector<double> Proba_Moves, Proba_Moves_Acc_Sum;  
   public:   

   ///
   move_set(polymorphic_random_generator & R): RNG(R) { Proba_Moves.push_back(0); }

   /** 
    * Add move M with its probability of being proposed.
    * NB : the PropositionProbability needs to be >0 but does not need to be 
    * normalized. Normalization is automatically done with all the added moves
    * before starting the run
    *
    * WARNING : the pointer is deleted automatically by the MC class at destruction. 
    */
   template <typename MoveType>
    void add (MoveType *M, std::string name, double PropositionProbability) {
     this->push_back( move<MCSignType> (M) );
     assert(PropositionProbability >=0);
     Proba_Moves.push_back(PropositionProbability);
     names_.push_back(name);
     normaliseProba();// ready to run after each add !
    }

   template <typename MoveType>
    void add (boost::shared_ptr<MoveType> sptr, std::string name, double PropositionProbability) {
     this->push_back( move<MCSignType> (sptr) );
     assert(PropositionProbability >=0);
     Proba_Moves.push_back(PropositionProbability);
     names_.push_back(name);
     normaliseProba();// ready to run after each add !
    }

   /**
    *  - Picks up one of the move at random (weighted by their proposition probability), 
    *  - Call Try method of that move
    *  - Returns the metropolis ratio R (see move concept). 
    *    The sign ratio returned by the try method of the move is kept.
    */
   double Try() {
    assert( Proba_Moves_Acc_Sum.size()>0);
    // Choice of move with its probability
    double proba = RNG(); assert(proba>=0);
    current_move_number =0; while (proba >= Proba_Moves_Acc_Sum[current_move_number] ) { current_move_number++;}
    assert(current_move_number>0); assert(current_move_number<=this->size());
    current_move_number--;
    current =  & (*this)[current_move_number];
#ifdef TRIQS_TOOLS_MC_DEBUG 
    std::cerr << "*******************************************************"<< std::endl;
    std::cerr << "Name of the proposed move: " << name_of_currently_selected() << std::endl;
    std::cerr <<"  Proposition probability = "<<proba<<std::endl;
#endif
    MCSignType rate_ratio = current->Try();
    if (!std::isfinite(std::abs(rate_ratio))) 
     TRIQS_RUNTIME_ERROR<<"Monte Carlo Error : the rate is not finite in move "<<name_of_currently_selected();
    double abs_rate_ratio = std::abs(rate_ratio);
#ifdef TRIQS_TOOLS_MC_DEBUG
    std::cerr << " Metropolis ratio " << rate_ratio<<". Abs(Metropolis ratio) " <<abs_rate_ratio << std::endl;
#endif
    assert ((abs_rate_ratio>=0)); 
    try_sign_ratio = ( abs_rate_ratio> 1.e-14 ? rate_ratio/abs_rate_ratio : 1); // keep the sign
    return abs_rate_ratio;
   }

   /**
    *  Accept the move previously selected and tried.
    *  Returns the Sign computed as, if M is the move : 
    *  Sign = sign (M.Try()) * M.Accept()
    */
   MCSignType Accept() { 
    MCSignType accept_sign_ratio =  current->Accept();
    // just make sure that accept_sign_ratio is a sign!
    assert(std::abs(std::abs(accept_sign_ratio)-1.0) < 1.e-10);
#ifdef TRIQS_TOOLS_MC_DEBUG
    std::cerr.setf(std::ios::scientific, std::ios::floatfield);
    std::cerr<<" ... Move accepted"<<std::endl;
    std::cerr<<"   try_sign_ratio = "  << try_sign_ratio   <<std::endl;
    std::cerr<<"   accept_sign_ratio = "<< accept_sign_ratio <<std::endl;
    std::cerr<<"   their product  =  "<< try_sign_ratio* accept_sign_ratio <<std::endl;
#endif
    return try_sign_ratio * accept_sign_ratio;
   }

   /**  Reject the move :  Call the Reject() method of the move previously selected
   */
   void Reject() { 
#ifdef TRIQS_TOOLS_MC_DEBUG
    std::cerr<<" ... Move rejected"<<std::endl;
#endif
    current->Reject();
   } 

   /// Pretty printing of the acceptance probability of the moves. 
   void print (triqs::utility::report_stream & report, mpi::communicator const & c, std::string name="", std::string decal="") const { 
    report <<decal <<"Move set : "<<name <<std::endl;
    for (unsigned int u =0; u< this->size(); ++u)
     (*this)[u].print(report,c,names_[u],decal+std::string("  ")); 
   }

   protected:
   MCSignType try_sign_ratio;

   void normaliseProba() { // Computes the normalised accumulated probability
    if (this->size() ==0)  TRIQS_RUNTIME_ERROR<<" no moves registered";
    double acc = 0; 
    Proba_Moves_Acc_Sum.clear();
    for (unsigned int u = 0; u<Proba_Moves.size(); ++u) acc+=Proba_Moves[u];assert(acc>0);
    for (unsigned int u = 0; u<Proba_Moves.size(); ++u) Proba_Moves_Acc_Sum.push_back(Proba_Moves[u]/acc);
    for (unsigned int u = 1; u<Proba_Moves_Acc_Sum.size(); ++u) Proba_Moves_Acc_Sum[u] += Proba_Moves_Acc_Sum[u-1];
    assert(std::abs(Proba_Moves_Acc_Sum[Proba_Moves_Acc_Sum.size()-1] -1)<1.e-13);
    Proba_Moves_Acc_Sum[Proba_Moves_Acc_Sum.size()-1] += 0.001; 
    // I shift the last proba acc so that even if random number in onecycle is 1 it is below that bound
    assert(Proba_Moves_Acc_Sum.size()==this->size()+1);
   }

   std::string name_of_currently_selected() const { return names_[current_move_number];} 
  };// class move_set

}}// end namespace
#endif

