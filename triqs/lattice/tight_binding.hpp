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
#ifndef TRIQS_LATTICE_TIGHTBINDINGS_H
#define TRIQS_LATTICE_TIGHTBINDINGS_H
#include "bravais_lattice_and_brillouin_zone.hpp"
#include <boost/unordered_map.hpp>

namespace triqs { namespace lattice_tools { 

 /**
   For tightbinding Hamiltonian with fully localised orbitals
   Model the ShortRangeFunctionOnBravaisLattice concept.
   Overlap between orbital is taken as unit matrix. 
  */
 class tight_binding { 
  public : 
   typedef std::vector<int>                                      arg_type;
   typedef matrix_view <dcomplex>                                   return_type;
   typedef matrix <dcomplex>                                        return_construct_type;
   typedef boost::unordered_map<arg_type, array<dcomplex,2> >   map_type;

   tight_binding(bravais_lattice const & bl__, map_type const & t_r);
   tight_binding(bravais_lattice const & bl__,boost::python::object dct);
   bravais_lattice const & lattice() const { return bl_;}
   size_t n_bands() const { return bl_.n_orbitals();}
   map_type::const_iterator begin() const { return tr.begin();} 
   map_type::const_iterator end()   const { return tr.end();} 

  protected:
   bravais_lattice bl_;
   map_type tr;
   void check();
 };

 /**
   Factorized version of hopping (for speed)
   k_in[:,n] is the nth vector
   In the result, R[:,:,n] is the corresponding hopping t(k)
  */
 array_view <dcomplex,3> hopping_stack (tight_binding const & TB, array_view<double,2> const & k_stack);
 // only the default array have a boost python converter....
 // not optimal ordering here . But why is this in FORTRAN ?
 //array_view <dcomplex,3> hopping_stack (tight_binding const & TB, array_view<double,2, arrays::TraversalOrderFortran> const & k_stack);

 std::pair<array<double,1>, array<double,2> > dos(tight_binding const & TB, size_t nkpts, size_t neps); 
 std::pair<array<double,1>, array<double,1> > dos_patch(tight_binding const & TB, const array<double,2> & triangles, size_t neps, size_t ndiv);
 array_view<double,2> energies_on_bz_path(tight_binding const & TB, K_view_type K1, K_view_type K2, size_t n_pts);
 array_view<double,2> energies_on_bz_grid(tight_binding const & TB, size_t n_pts);

}}

#endif


