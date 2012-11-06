
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

#ifndef TRIQS_TOOLS_MC_MEASURE2_H
#define TRIQS_TOOLS_MC_MEASURE2_H 

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/function.hpp>
#include <boost/mpi.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/concept_check.hpp>
#include <map>
#include <triqs/utility/exceptions.hpp>

namespace triqs { namespace mc_tools { 
 namespace mpi=boost::mpi;
 namespace BLL = boost::lambda;
 
 template <class X, typename MCSignType> struct IsMeasure {
  BOOST_CONCEPT_USAGE(IsMeasure)
  {
   i.accumulate(MCSignType(1.0));    
   i.collect_results(*c);
  }
  private:
  X i; boost::mpi::communicator const * c;
 };

 //--------------------------------------------------------------------

 template<typename MCSignType>  
  class mcmeasure  { 
   boost::shared_ptr< void > impl_;
   boost::function<void (MCSignType const & ) > accumulate_;
   uint64_t count_;

   public :
   boost::function<void (boost::mpi::communicator const & )> collect_results;

   mcmeasure():impl_(), count_(0) {} // needed to make vector of these...

   template<typename MeasureType> mcmeasure ( MeasureType * p) : 
    impl_(p),
    accumulate_(BLL::bind(&MeasureType::accumulate,p,BLL::_1)),
    collect_results( BLL::bind(&MeasureType::collect_results, p,BLL::_1)),
    count_(0)
   { BOOST_CONCEPT_ASSERT((IsMeasure<MeasureType,MCSignType>)); }

   template<typename MeasureType> mcmeasure ( boost::shared_ptr<MeasureType> sptr) :
    impl_(sptr),
    accumulate_(BLL::bind(&MeasureType::accumulate,sptr.get(),BLL::_1)),
    collect_results(BLL::bind(&MeasureType::collect_results,sptr.get(),BLL::_1)),
    count_(0)
   { BOOST_CONCEPT_ASSERT((IsMeasure<MeasureType,MCSignType>)); }

   void accumulate(MCSignType signe){ assert(impl_); count_++; accumulate_(signe); }
 
   uint64_t count() const { return count_;}
};

//--------------------------------------------------------------------

template<typename MCSignType>
class measure_set : public std::map<std::string, mcmeasure<MCSignType> > {
 typedef std::map<std::string, mcmeasure<MCSignType> > BaseType;
 typedef mcmeasure<MCSignType> measure_type;
 public : 
 typedef typename BaseType::iterator iterator;
 typedef typename BaseType::const_iterator const_iterator;

 /**
  * Register the Measure M with a name
  * WARNING : the pointer is deleted automatically by the class at destruction. 
  */
 template<typename T>
  void insert (T *M, std::string const & name) {
   if (has(name)) TRIQS_RUNTIME_ERROR <<"measure_set : insert : measure '"<<name<<"' already inserted";
   BaseType::insert(std::make_pair(name, measure_type (M) )); 
  } 

 template<typename T>
  void insert (boost::shared_ptr<T> sptr, std::string const & name) {
   if (has(name)) TRIQS_RUNTIME_ERROR <<"measure_set : insert : measure '"<<name<<"' already inserted";
   BaseType::insert(std::make_pair(name, measure_type (sptr) ));
  }

 bool has(std::string const & name) const { return BaseType::find(name) != BaseType::end(); }

 measure_type & operator[](std::string const & name) {
  if (!has(name)) throw std::out_of_range("No result found with the name: " + name);
  return BaseType::find(name)->second;
 }

 measure_type const & operator[](std::string const & name) const {
  if (!has(name)) throw std::out_of_range("No result found with the name: " + name);
  return BaseType::find(name)->second;
 }

 ///
 void accumulate( MCSignType & signe) { for (iterator it= this->begin(); it != this->end(); ++it) it->second.accumulate(signe); }

 ///
 std::vector<std::string> names() const { 
  std::vector<std::string> res; 
  for (const_iterator it= this->begin(); it !=this->end(); ++it) res.push_back(it->first); 
  return res;
 }

 // gather result for all measure, on communicator c
 void collect_results (boost::mpi::communicator const & c ) {
  for (typename BaseType::iterator it = this->begin(); it != this->end(); ++it) it->second.collect_results(c);
 }

};

}}// end namespace
#endif

