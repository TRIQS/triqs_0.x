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
#ifndef TRIQS_TOOLS_MC_MEASURE2_H
#define TRIQS_TOOLS_MC_MEASURE2_H

#include <functional>
#include <boost/mpi.hpp>
#include <map>
#include <triqs/utility/exceptions.hpp>

namespace triqs { namespace mc_tools {

 // mini concept checking
 template<typename MCSignType, typename T, typename Enable=void> struct has_accumulate : std::false_type {};
 template<typename MCSignType, typename T> struct has_accumulate <MCSignType, T, decltype(std::declval<T>().accumulate(MCSignType()))> : std::true_type {};

 template<typename T, typename Enable=void> struct has_collect_result : std::false_type {};
 template<typename T> struct has_collect_result < T, decltype(std::declval<T>().collect_results(std::declval<boost::mpi::communicator>()))> : std::true_type {};

 //--------------------------------------------------------------------

 template<typename MCSignType>
  class mcmeasure {

   std::shared_ptr<void> impl_;
   std::function<mcmeasure()> clone_;
   size_t hash_;
   std::string type_name_;

   std::function<void (MCSignType const & ) > accumulate_;
   std::function<void (boost::mpi::communicator const & )> collect_results_;
   uint64_t count_;

   template<typename MeasureType>
    void deleg (MeasureType * p) {
    impl_= std::shared_ptr<void> (p);
    accumulate_ = [p](MCSignType const & x) { p->accumulate(x);};
    count_ = 0;
    hash_ = typeid(MeasureType).hash_code();
    type_name_ =  typeid(MeasureType).name();
    collect_results_ = [p] ( boost::mpi::communicator const & c) { p->collect_results(c);};
    clone_ =  [p]() { return MeasureType(*p);} ;
    //h5_w       = [obj](h5::group F, std::string const &Name)->void { h5_write(F,Name, *obj);};
   }

   template<typename MeasureType> void check_code() const {
    if (typeid(MeasureType).hash_code() != hash_)
     TRIQS_RUNTIME_ERROR << "Trying to retrieve a measure of type "<< typeid(MeasureType).name() << " from a measure of type "<< type_name_;
   };

   public :

   template<typename MeasureType>
    mcmeasure (MeasureType && p) {
     static_assert( has_accumulate<MCSignType,MeasureType>::value, " Measure has no accumulate method !");
     static_assert( has_collect_result<MeasureType>::value, " Measure has no collect_results method !");
     deleg( new typename std::remove_reference<MeasureType>::type(std::forward<MeasureType>(p)));
    }

   // Value semantics. Everyone at the end call move = ...
   mcmeasure(mcmeasure const &rhs) {*this = rhs;}
   mcmeasure(mcmeasure &rhs) {*this = rhs;} // or it will use the template  = bug
   mcmeasure(mcmeasure && rhs) { *this = std::move(rhs);}
   mcmeasure & operator = (mcmeasure const & rhs) { *this = rhs.clone_(); return *this;}
   mcmeasure & operator = (mcmeasure && rhs) {
    using std::swap;
    swap(impl_,rhs.impl_); swap(accumulate_, rhs.accumulate_); swap(count_, rhs.count_); swap( hash_, rhs.hash_); swap(type_name_,rhs.type_name_);
    swap(collect_results_,rhs.collect_results_); swap(clone_,rhs.clone_);
    return *this;
   }

   void accumulate(MCSignType signe){ assert(impl_); count_++; accumulate_(signe); }
   void collect_results (boost::mpi::communicator const & c ) { collect_results_(c);}

   uint64_t count() const { return count_;}
   size_t hash_code() const { return hash_;}
   template<typename MeasureType> MeasureType       & get()       { check_code<MeasureType>(); return *(static_cast<MeasureType *>(impl_.get())); }
   template<typename MeasureType> MeasureType const & get() const { check_code<MeasureType>(); return *(static_cast<MeasureType const *>(impl_.get())); }
  };

 //--------------------------------------------------------------------

 template<typename MCSignType>
  class measure_set  {
   typedef mcmeasure<MCSignType> measure_type;
   std::map<std::string, mcmeasure<MCSignType>> m_map;
   public :

   /**
    * Register the Measure M with a name
    */
   template<typename MeasureType>
    void insert (MeasureType && M, std::string const & name) {
     if (has(name)) TRIQS_RUNTIME_ERROR <<"measure_set : insert : measure '"<<name<<"' already inserted";
     m_map.insert(std::make_pair(name, measure_type (std::forward<MeasureType>(M))));
     // not implemented on gcc 4.6's stdlibc++ ?
     // m_map.emplace(name, measure_type (std::forward<MeasureType>(M)));
    }

   bool has(std::string const & name) const { return m_map.find(name) != m_map.end(); }

   void accumulate(MCSignType & signe) { for (auto & nmp : m_map) nmp.second.accumulate(signe); }

   std::vector<std::string> names() const {
    std::vector<std::string> res;
    for (auto & nmp : m_map) res.push_back(nmp.first);
    return res;
   }

   // gather result for all measure, on communicator c
   void collect_results (boost::mpi::communicator const & c ) { for (auto & nmp : m_map) nmp.second.collect_results(c); }

   // access to the measure, given its type, with dynamical type check
   template<typename MeasureType>
    MeasureType & get_measure(std::string const & name) {
     auto it = m_map.find (name);
     if (it == m_map.end()) TRIQS_RUNTIME_ERROR << " Measure " << name << " unknown";
     return it->template get<MeasureType>();
    }

   template<typename MeasureType>
    MeasureType const & get_measure(std::string const & name) const {
     auto it = m_map.find (name);
     if (it == m_map.end()) TRIQS_RUNTIME_ERROR << " Measure " << name << " unknown";
     return it->template get<MeasureType>();
    }
  };

}}// end namespace
#endif

