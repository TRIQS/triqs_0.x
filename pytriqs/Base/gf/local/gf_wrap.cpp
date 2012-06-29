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
#include<Python.h>
#include <triqs/gf/local/gf.hpp>
#include <boost/python/def.hpp>
#include <triqs/arrays/python/converters.hpp>

using namespace boost::python;
//namespace tqa=triqs::arrays;

using namespace triqs::gf;
using local::tail_view;

BOOST_PYTHON_MODULE(_pytriqs_GF2) {

 int _r = _import_array();assert(_r==0); // to be refined...
// tqa::register_boost_converters();

 //register_exception_translator<std::string>(translatorString);
 //triqs::python_tools::register_converter< triqs::arrays::array<double,2> > ();
 triqs::python_tools::register_converter< triqs::arrays::array_view<std::complex<double>,3> > ();

 docstring_options doc_options;
 doc_options.disable_cpp_signatures();
 doc_options.disable_py_signatures();

 // **********   ******************

 enum_<statistic_enum>("GF_Statistic") .value("Boson", Boson) .value("Fermion",Fermion) ;

 // **********   ******************

 class_<domains::matsubara_freq>("DomainMatsubaraFrequency", init<double, statistic_enum>())
  .add_property("Beta",&domains::matsubara_freq::beta)
  .add_property("Statistic",&domains::matsubara_freq::statistic)
  ;

 class_<meshes::matsubara_freq> ("MeshMatsubaraFrequency", init<double, statistic_enum, size_t>())
  .def("domain",&meshes::matsubara_freq::domain,return_internal_reference<1>()) // TO BE TESTED !! 
  .def("size",&meshes::matsubara_freq::size)
  ;

 // **********  tail  ******************

 class_<tail_view>("TailGF", init<tail_view::data_type,long>() )
  ;

};
