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

#include "bravais_lattice_and_brillouin_zone.hpp"
#include "tight_binding.hpp"
#include <boost/python.hpp>
#include <boost/python/docstring_options.hpp>
#include <triqs/utility/compiler_details.hpp>
#include <triqs/arrays/python/converters.hpp>
#include <triqs/python_tools/converters/pair.hpp> 
#include <triqs/python_tools/converters/unordered_map.hpp> 

using namespace boost::python;

namespace triqs { namespace lattice_tools {  

 BOOST_PYTHON_MODULE(_pytriqs_LatticeTools) {

int _r = _import_array();assert(_r==0);
  tqa::register_boost_converters();
  triqs::python_tools::register_converter< std::pair<array<double,1>, tqa::array<double,1> > >();
  triqs::python_tools::register_converter< std::pair<array<double,1>, tqa::array<double,2> > >();

  //triqs::python_tools::register_converter< tqa::array_view<double,2> >();
  //triqs::python_tools::register_converter< tqa::array_view<double,2, tqa::Option::Fortran> >();
  //triqs::python_tools::register_converter< tqa::vector_view<double> >();

  docstring_options doc_options;
  doc_options.disable_py_signatures();
  doc_options.disable_cpp_signatures();

  char bravais_lattice_doc_init[] = 
   " :param Units: a list of 3d vectors representing the units of the lattice in R^3 \n"
   "               the dimension of the lattice is the number of these vectors.\n"
   " :param Orbital_Positions: Position of the orbitals close to the origin. Default = (0,0,0)";

  class_<bravais_lattice>("bravais_lattice", 
    init<object , object  >
    //init<array_view<double,2>, object  >
    ( bravais_lattice_doc_init,
      ( boost::python::arg("Units"), 
	boost::python::arg("Orbital_Positions")=boost::python::list(make_tuple(make_tuple(0.0,0.0,0.0)))
      )
    ))
   .def("lattice_to_real_coordinates",&bravais_lattice::lattice_to_real_coordinates,
     "Transform into real coordinates.\n"
     ":param x_in: coordinates in the basis of the unit vectors\n"
     ":rtype: Output coordinates in R^3 (real coordinates) as an array")
   .def_readonly("dim",  &bravais_lattice::dim, "Dimension of the lattice")
   .def_readonly("n_orbitals",&bravais_lattice::n_orbitals, "Number of orbitals")
   ;

  char tight_binding_doc_init[] = 
   " :param BravaisLattice: The underlying Bravais lattice\n"
   " :param Hopping: A dictionnary coding the tight-binding hopping\n\n"
   "        * keys are displacement R on the lattice : d-dimensional vectors of integers\n"
   "        * values are matrices t_ab(R), as 2d numpy or list of list (anything from which numpy can make a 2d array)\n";  

  class_<tight_binding>("tight_binding", 
    init <bravais_lattice, object >(tight_binding_doc_init, ( arg("BravaisLattice"), arg("Hopping") ))
    )
   ;

  def("dos_patch",dos_patch);
  def("dos", dos);  
  def ("energies_on_bz_grid",energies_on_bz_grid);
  def ("energies_on_bz_path",energies_on_bz_path);
  def ("hopping_stack", hopping_stack);

 };

}}


