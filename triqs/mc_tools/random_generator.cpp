
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

#include "random_generator.hpp"
#include "./MersenneRNG.hpp"
#include <boost/random.hpp>
//#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
//#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/ranlux.hpp>
#include <sstream>
#include <boost/preprocessor/seq.hpp>
#include <boost/preprocessor/control/if.hpp>
#define AS_STRING(X) AS_STRING2(X)
#define AS_STRING2(X) #X

// List of All available Boost random number generator
#define RNG_LIST (mt19937)(mt11213b)\
(lagged_fibonacci607) (lagged_fibonacci1279) (lagged_fibonacci2281) (lagged_fibonacci3217) (lagged_fibonacci4423)\
(lagged_fibonacci9689) (lagged_fibonacci19937) (lagged_fibonacci23209) (lagged_fibonacci44497) (ranlux3)


namespace triqs { 
 namespace mc_tools { 

 typedef boost::generator <double> gen_type;

 inline gen_type * choose_gen(std::string const & RandomGeneratorName, std::size_t seed_, boost::shared_ptr<void> &ptr) { 

  if (RandomGeneratorName=="") return new gen_type( mc_tools::RandomGenerators::RandMT (seed_));

  // now boost random number generators
#define DRNG(r,data,XX) if (RandomGeneratorName==AS_STRING(XX)) { \
 boost::shared_ptr<boost::XX> localptr = boost::make_shared<boost::XX>(seed_);\
 ptr = localptr; boost::uniform_real<> dis;\
 return new gen_type( boost::variate_generator<boost::XX&, boost::uniform_real<> >(*localptr,dis));}

  BOOST_PP_SEQ_FOR_EACH(DRNG,~,RNG_LIST)

   TRIQS_RUNTIME_ERROR<<"The random generator "<<RandomGeneratorName<<" is not recognized";
 }

 //---------------------------------------------

 random_generator::random_generator(std::string const & RandomGeneratorName, std::size_t seed_ ) : 
  rng_ptr(), gen (choose_gen(RandomGeneratorName,seed_,rng_ptr)), name(RandomGeneratorName),seed(seed_) {}

 //---------------------------------------------

 random_generator::random_generator( random_generator const & p) :
  rng_ptr(), gen (choose_gen(p.name,p.seed,rng_ptr)), name (p.name),seed(p.seed) {}

 //---------------------------------------------

 std::string random_generator::random_generator_names(std::string const & sep) { 
#define PR(r,sep,p,XX) BOOST_PP_IF(p,+ sep +,) std::string(AS_STRING(XX))   
  return BOOST_PP_SEQ_FOR_EACH_I (PR,sep,RNG_LIST);  
 }
}
}
