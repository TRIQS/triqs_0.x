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
#include <triqs/arrays/array.hpp>
#include <triqs/arrays/matrix.hpp>
#include <triqs/arrays/vector.hpp>
#include <triqs/utility/exceptions.hpp>
#include <map>
#include <string>
#include <boost/unordered_map.hpp>

namespace triqs { namespace lattice {

 namespace tqa = triqs::arrays;
 using tqa::array; using tqa::array_view; using tqa::matrix_view;using tqa::matrix;
 typedef std::complex<double> dcomplex;

 /**
  */
 class bravais_lattice {
  public :
   typedef tqa::vector <double> point_t;
   typedef tqa::matrix<double> units_t;
   typedef boost::unordered_map<std::string, point_t> orbital_t; // name -> position in space

   bravais_lattice( units_t const & units__, orbital_t const & orbitals__);
   size_t n_orbitals() const {return orbitals_.size();}
   units_t const & units() const { return units_;}
   size_t dim() const { return dim_; }

   /***
    * Transform into real coordinates.
    * @param[in] x : coordinates in the basis :unit_i
    * @return  Coordinates in R^3 ("real" coordinates)
    */
   point_t lattice_to_real_coordinates(point_t const & x) const;

   protected :
   units_t units_;
   size_t dim_;
   orbital_t orbitals_;
   void cons_deleg(tqa::array<double,2> const & units__);
 };

 /**
  * Brillouin Zone class ...
  */
 class brillouin_zone {
  bravais_lattice lattice_;
  tqa::matrix<double> K_reciprocal, K_reciprocal_inv;
  public :

  typedef tqa::vector <double> point_t;
  brillouin_zone( bravais_lattice const & bl_);
  bravais_lattice lattice() const { return lattice_;}

  /**
   * Transform into real coordinates.
   * @param[in] k : coordinates in the basis (K_reciprocal_i)
   * @return  Coordinates in R^3 ("real" coordinates)
   */
  point_t lattice_to_real_coordinates (point_t const & k) const;

  /// Inverse of latt_to_real_k
  point_t real_to_lattice_coordinates (point_t const & k) const;
 };

}}
#endif
