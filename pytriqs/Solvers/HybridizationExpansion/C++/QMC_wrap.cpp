
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

#include <triqs/utility/h5_exceptions.hpp>
#include "MC.hpp"
#include <boost/python/return_internal_reference.hpp>
#include <GreenFunctions/GF_Bloc_ImTime.hpp>
#include "Hloc.hpp"

using namespace boost::python;

namespace MC_Hybridization_Matsu {
 void solve (boost::python::object );
};

BOOST_PYTHON_MODULE(_pytriqs_Solver_HybridizationExpansion) {

 triqs::utility::register_h5_exception();

 docstring_options doc_options;
 doc_options.disable_py_signatures();

 class_<Hloc>("Hloc",init<int,int,python::dict,python::dict,python::list,python::object,int>())
  .def ("__repr__",&Hloc::print)
  ;

 def ("MC_solve",&MC_Hybridization_Matsu::solve);

 def ("Random_Generators_Available", &triqs::mc_tools::polymorphic_random_generator::random_generator_names);

};
