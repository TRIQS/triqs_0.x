/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012 by M. Ferrero, O. Parcollet
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
#ifndef TRIQS_GF_DOMAIN_MATSUBARA_H
#define TRIQS_GF_DOMAIN_MATSUBARA_H
#include "../tools.hpp"

namespace triqs { namespace gf { 

 /// The domain
 template<bool IsComplex>
 struct matsubara_domain {
  typedef typename mpl::if_c<IsComplex, std::complex<double>, double>::type point_t;
  double beta;
  statistic_enum statistic;
  matsubara_domain (double Beta=1, statistic_enum s = Fermion): beta(Beta), statistic(s){}
  bool operator == (matsubara_domain const & D) const { return ((std::abs(beta - D.beta)<1.e-15) && (statistic == D.statistic));}

  /// Write into HDF5
  friend void h5_write (tqa::h5::group_or_file fg, std::string subgroup_name, matsubara_domain const & d) {
   tqa::h5::group_or_file gr =  fg.create_group(subgroup_name);
   h5_write(gr,"Beta",d.beta);
   h5_write(gr,"Statistic",(d.statistic==Fermion ? 'F' : 'B'));
  }

  /// Read from HDF5
  friend void h5_read  (tqa::h5::group_or_file fg, std::string subgroup_name, matsubara_domain & d){
   tqa::h5::group_or_file gr = fg.open_group(subgroup_name);
   double beta; char statistic;
   h5_read(gr,"Beta",beta);
   h5_read(gr,"Statistic",statistic);
   d = matsubara_domain(beta,(statistic=='F' ? Fermion : Boson));
  }
 
  //  BOOST Serialization
  friend class boost::serialization::access;
  template<class Archive>
   void serialize(Archive & ar, const unsigned int version) {
    ar & boost::serialization::make_nvp("Beta",beta);
    ar & boost::serialization::make_nvp("Statistic",statistic);
   }

 };

}}
#endif

