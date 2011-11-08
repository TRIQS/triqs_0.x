
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

#ifndef TRIQS_LATTICE_BRAVAIS_LATTICE_H
#define TRIQS_LATTICE_BRAVAIS_LATTICE_H
#include <boost/python.hpp>
#include <boost/unordered_map.hpp>
#include <triqs/arrays/array.hpp>
#include <triqs/arrays/matrix.hpp>
#include <triqs/arrays/vector.hpp>
#include <triqs/utility/exceptions.hpp>
#include <map>
#include <string>

namespace triqs { namespace lattice_tools { 

 namespace tqa = triqs::arrays; 
 using tqa::array; using tqa::array_view; using tqa::matrix_view;using tqa::matrix;

 typedef tqa::vector <double> R_type;
 typedef tqa::vector <double> K_type;
 typedef tqa::vector_view<double> R_view_type;
 typedef tqa::vector_view<double> K_view_type;
 typedef std::complex<double> dcomplex;

 class bravais_lattice { 
  public : 
   typedef tqa::matrix<double> units_type;
   typedef boost::unordered_map<std::string, R_type> orbital_type; // name -> position in space
   
   bravais_lattice( units_type const & units__, orbital_type const & orbitals__); 
   bravais_lattice( tqa::array<double,2> const & units__, boost::python::object orbitals__); 
   size_t n_orbitals() const {return orbitals_.size();}
   units_type const & units() const { return units_;}
   size_t dim() const { return dim_; }   
   /***
    * Transform into real coordinates.
    * @param[in] x : coordinates in the basis :unit_i
    * @return  Coordinates in R^3 ("real" coordinates) 
    */
   R_view_type lattice_to_real_coordinates(R_type const & x) const;

   protected :
   units_type units_;
   size_t dim_;
   orbital_type orbitals_;
   void cons_deleg(tqa::array<double,2> const & units__);
 };

 class brillouin_zone { 
  bravais_lattice lattice_;
  tqa::matrix<double> K_reciprocal, K_reciprocal_inv;
  public : 

  brillouin_zone( bravais_lattice const & bl_);
  bravais_lattice lattice() const { return lattice_;}

  /***
   * Transform into real coordinates.
   * @param[in] k : coordinates in the basis (K_reciprocal_i)
   * @return  Coordinates in R^3 ("real" coordinates) 
   */
  K_view_type lattice_to_real_coordinates (K_view_type const & k) const; 

  /// Inverse of latt_to_real_k
  K_view_type real_to_lattice_coordinates (K_view_type const & k) const;
 };

}}
#endif
