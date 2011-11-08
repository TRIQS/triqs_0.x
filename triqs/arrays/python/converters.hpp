
/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011 by O. Parcollet
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

#ifndef TRIQS_ARRAYS_BOOST_CONVERTER_H
#define TRIQS_ARRAYS_BOOST_CONVERTER_H

#ifndef TRIQS_ARRAYS_WITH_PYTHON_SUPPORT
#error "You must define the macro TRIQS_ARRAYS_WITH_PYTHON_SUPPORT to use Python interface"
#endif

//#include <boost/python/slice.hpp>
#include "../python/numpy_interface.hpp"
#include <complex>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>

/* ***********************************************

   Mainly for debug : in real code, it is better to list explicitely the converters needed.
   Function to register the convert between numpy object and triqs_arrays into boost.python
Usage : (in the main, before launching the interpreter)
#include <triqs/arrays/python/converters.hpp>"
triqs::arrays::register_boost_converters();

 ****************************************************/

// List all the types for which a converter will be generated.
#define TYPES_FOR_WHICH_CONVERTER_IS_NEEDED (long)(std::complex<double>)(double) 

namespace triqs { namespace arrays { 

 template<typename T>
  inline void register_boost_converters_for_one_type() {
#define ONE(z,N,unused)\
   triqs::python_tools::register_converter< triqs::arrays::array<T,BOOST_PP_INC(N),triqs::arrays::Option::C> > ();\
   triqs::python_tools::register_converter< triqs::arrays::array<T,BOOST_PP_INC(N),triqs::arrays::Option::Fortran> > ();\
   triqs::python_tools::register_converter< triqs::arrays::array_view<T,BOOST_PP_INC(N),triqs::arrays::Option::C> > ();\
   triqs::python_tools::register_converter< triqs::arrays::array_view<T,BOOST_PP_INC(N),triqs::arrays::Option::Fortran> > ();
   BOOST_PP_REPEAT(ARRAY_NRANK_MAX, ONE, nil);
#undef ONE

#define CONVERT_ALSO_MATRIX_VECTOR 
#ifdef CONVERT_ALSO_MATRIX_VECTOR
   triqs::python_tools::register_converter< triqs::arrays::matrix<T,triqs::arrays::Option::C> > ();
   triqs::python_tools::register_converter< triqs::arrays::matrix<T,triqs::arrays::Option::Fortran> > ();
   triqs::python_tools::register_converter< triqs::arrays::matrix_view<T,triqs::arrays::Option::C> > ();
   triqs::python_tools::register_converter< triqs::arrays::matrix_view<T,triqs::arrays::Option::Fortran> > ();
   triqs::python_tools::register_converter< triqs::arrays::vector<T,triqs::arrays::Option::C> > ();
   triqs::python_tools::register_converter< triqs::arrays::vector_view<T,triqs::arrays::Option::C> > ();
#endif
  }

 inline void register_boost_converters() {
  triqs::python_tools::register_converter< triqs::arrays::range> ();
#define ONETYPE(r, data, elem) register_boost_converters_for_one_type<elem>();
  BOOST_PP_SEQ_FOR_EACH(ONETYPE, nil , TYPES_FOR_WHICH_CONVERTER_IS_NEEDED); 
#undef ONETYPE
 }

}}//namespace triqs::arrays 

#undef TYPES_FOR_WHICH_CONVERTER_IS_NEEDED
#endif 

