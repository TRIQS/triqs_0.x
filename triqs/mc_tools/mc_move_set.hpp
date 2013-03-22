/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011-2013 by M. Ferrero, O. Parcollet
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
#include <triqs/utility/report_stream.hpp>
#include <triqs/utility/exceptions.hpp>
#include <functional>
#include <boost/mpi.hpp>
#include "./random_generator.hpp"

namespace triqs { namespace mc_tools { 

 // mini concept checking
 template<typename MCSignType, typename T, typename Enable=MCSignType> struct has_attempt : std::false_type {};
 template<typename MCSignType, typename T> struct has_attempt <MCSignType, T, decltype(MCSignType(std::declval<T>().attempt()))> : std::true_type {};

 template<typename MCSignType, typename T, typename Enable=MCSignType> struct has_accept : std::false_type {};
 template<typename MCSignType, typename T> struct has_accept <MCSignType,T, decltype(std::declval<T>().accept())> : std::true_type {};

 template<typename T, typename Enable=void> struct has_reject : std::false_type {};
 template<typename T> struct has_reject<T, decltype(std::declval<T>().reject())> : std::true_type {};

 //----------------------------------
 
 template<typename MCSignType> class move_set;

 template<typename MCSignType> 
  class move {

   std::shared_ptr<void> impl_;
   std::function<move()> clone_; 
   size_t hash_; 
   std::string type_name_;

   std::function<MCSignType()> attempt_, accept_;
   std::function<void()> reject_;
   uint64_t NProposed, Naccepted;
   double acceptance_rate_;
   move_set<MCSignType> * mset_ptr_;

   static move_set<MCSignType> * _mk_ptr( move_set<MCSignType> * x) { return x;}
   template<class T> static move_set<MCSignType> * _mk_ptr( T * x) { return NULL;}

   template<typename MoveType> 
    void deleg (MoveType * p) {
     impl_= std::shared_ptr<void> (p);
     hash_ = typeid(MoveType).hash_code();
     type_name_ =  typeid(MoveType).name();
     clone_   = [p]() { return MoveType(*p);};
     attempt_ = [p]() { return p->attempt();};
     accept_  = [p]() { return p->accept();};
     reject_  = [p]() { p->reject();};
     NProposed=0;
     Naccepted=0;
     acceptance_rate_ =-1;
     //h5_w       = [obj](h5::group F, std::string const &Name)->void { h5_write(F,Name, *obj);};
     mset_ptr_ =  _mk_ptr(p);
    }

   template<typename MoveType> void check_code() const { 
    if (typeid(MoveType).hash_code() != hash_) 
     TRIQS_RUNTIME_ERROR << "Trying to retrieve a measure of type "<< typeid(MoveType).name() << " from a measure of type "<< type_name_;
   };

   public :

   template<typename MoveType> 
    move (MoveType && p ) { 
     static_assert( has_attempt<MCSignType,MoveType>::value, " Move has no attempt method !"); 
     static_assert( has_accept<MCSignType,MoveType>::value, " Move has no accept method !"); 
     static_assert( has_reject<MoveType>::value, " Move has no reject method !"); 
     deleg( new typename std::remove_reference<MoveType>::type(std::forward<MoveType>(p)));
    }

   // Value semantics. Everyone at the end call move = ...
   move(move const &rhs) {*this = rhs;}  
   move(move &rhs) {*this = rhs;} // to avoid clash with tempalte construction 
   move(move && rhs) { *this = std::move(rhs);}
   move & operator = (move const & rhs) { *this = rhs.clone_(); return *this;}
   move & operator = (move && rhs) { 
    using std::swap; 
    swap(impl_,rhs.impl_); swap( hash_, rhs.hash_); swap(type_name_,rhs.type_name_); swap(clone_,rhs.clone_);
    swap(attempt_,rhs.attempt_); swap(accept_, rhs.accept_); swap(reject_, rhs.reject_);
    swap(NProposed,rhs.NProposed);swap(Naccepted,rhs.Naccepted); swap(acceptance_rate_,rhs.acceptance_rate_);
    swap(mset_ptr_,rhs.mset_ptr_);
    return *this;
   }

   MCSignType attempt(){ NProposed++; return attempt_();}
   MCSignType accept() { Naccepted++; return accept_(); }
   void reject() { reject_(); }
   double acceptance_rate() const { return acceptance_rate_;}
   move_set<MCSignType> * mset_ptr() { return mset_ptr_;}

   void collect_statistics(boost::mpi::communicator const & c) {
    uint64_t nacc_tot=0, nprop_tot=1;
    boost::mpi::reduce(c, Naccepted, nacc_tot,  std::plus<uint64_t>(), 0);
    boost::mpi::reduce(c, NProposed, nprop_tot, std::plus<uint64_t>(), 0);
    acceptance_rate_ = nacc_tot/static_cast<double>(nprop_tot);
   }

  };

 //--------------------------------------------------------------------

 /// A vector of (moves, proposition_probability), which is also a move itself
 template<typename MCSignType>
  class move_set  { 
   std::vector<move<MCSignType> > move_vec; 
   std::vector<std::string> names_;
   move<MCSignType> * current;
   size_t current_move_number;
   random_generator * RNG;
   std::vector<double> Proba_Moves, Proba_Moves_Acc_Sum;  
   public:   

   ///
   move_set(random_generator & R): RNG(&R) { Proba_Moves.push_back(0); }

   /** 
    * Add move M with its probability of being proposed.
    * NB : the proposition_probability needs to be >0 but does not need to be 
    * normalized. Normalization is automatically done with all the added moves
    * before starting the run
    */
   template <typename MoveType>
    void add (MoveType && M, std::string name, double proposition_probability) {
     move_vec.push_back( move<MCSignType> (std::forward<MoveType>(M)) );
     assert(proposition_probability >=0);
     Proba_Moves.push_back(proposition_probability);
     names_.push_back(name);
     normaliseProba();// ready to run after each add !
    }

   /**
    *  - Picks up one of the move at random (weighted by their proposition probability), 
    *  - Call attempt method of that move
    *  - Returns the metropolis ratio R (see move concept). 
    *    The sign ratio returned by the try method of the move is kept.
    */
   double attempt() {
    assert( Proba_Moves_Acc_Sum.size()>0);
    // Choice of move with its probability
    double proba = (*RNG)(); assert(proba>=0);
    current_move_number =0; while (proba >= Proba_Moves_Acc_Sum[current_move_number] ) { current_move_number++;}
    assert(current_move_number>0); assert(current_move_number<=move_vec.size());
    current_move_number--;
    current =  & move_vec[current_move_number];
#ifdef TRIQS_TOOLS_MC_DEBUG 
    std::cerr << "*******************************************************"<< std::endl;
    std::cerr << "Name of the proposed move: " << name_of_currently_selected() << std::endl;
    std::cerr <<"  Proposition probability = "<<proba<<std::endl;
#endif
    MCSignType rate_ratio = current->attempt();
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
    *  accept the move previously selected and tried.
    *  Returns the Sign computed as, if M is the move : 
    *  Sign = sign (M.attempt()) * M.accept()
    */
   MCSignType accept() { 
    MCSignType accept_sign_ratio =  current->accept();
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

   /**  reject the move :  Call the reject() method of the move previously selected
   */
   void reject() { 
#ifdef TRIQS_TOOLS_MC_DEBUG
    std::cerr<<" ... Move rejected"<<std::endl;
#endif
    current->reject();
   } 

   /// Pretty printing of the acceptance probability of the moves. 
   std::string get_statistics(boost::mpi::communicator const & c, int shift = 0) {
    std::ostringstream s;
    for (unsigned int u =0; u< move_vec.size(); ++u) {
     move_vec[u].collect_statistics(c);
     for(int i=0; i<shift; i++) s << " ";
     if (move_vec[u].mset_ptr()) {
      s << "Move set " << names_[u] << ": " << move_vec[u].acceptance_rate() << "\n";
      s << move_vec[u].mset_ptr()->get_statistics(c,shift+2);
     } else {
      s << "Move " << names_[u] << ": " << move_vec[u].acceptance_rate() << "\n";
     }
    }
    return s.str();
   }

   protected:
   MCSignType try_sign_ratio;

   void normaliseProba() { // Computes the normalised accumulated probability
    if (move_vec.size() ==0)  TRIQS_RUNTIME_ERROR<<" no moves registered";
    double acc = 0; 
    Proba_Moves_Acc_Sum.clear();
    for (unsigned int u = 0; u<Proba_Moves.size(); ++u) acc+=Proba_Moves[u];assert(acc>0);
    for (unsigned int u = 0; u<Proba_Moves.size(); ++u) Proba_Moves_Acc_Sum.push_back(Proba_Moves[u]/acc);
    for (unsigned int u = 1; u<Proba_Moves_Acc_Sum.size(); ++u) Proba_Moves_Acc_Sum[u] += Proba_Moves_Acc_Sum[u-1];
    assert(std::abs(Proba_Moves_Acc_Sum[Proba_Moves_Acc_Sum.size()-1] -1)<1.e-13);
    Proba_Moves_Acc_Sum[Proba_Moves_Acc_Sum.size()-1] += 0.001; 
    // I shift the last proba acc so that even if random number in onecycle is 1 it is below that bound
    assert(Proba_Moves_Acc_Sum.size()==move_vec.size()+1);
   }

   std::string name_of_currently_selected() const { return names_[current_move_number];} 
  };// class move_set

}}// end namespace
#endif

