
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

#ifndef PYARRAY_CONVERTER_H
#define PYARRAY_CONVERTER_H

#include "blitz_op.hpp"
#include "PyArray.hpp"

namespace triqs { 
namespace blitzext { 

 template<typename T>
  inline void myreg() {
   register_converter_python_array<T,1>();
   register_converter_python_array<T,2>();
   register_converter_python_array<T,3>();
   register_converter_python_array<T,4>();
   register_converter_python_array<T,5>();
   register_converter_python_array<T,6>();
  }

 inline void init_boost_array_converters() { 
  static bool initialized = false;
  if (initialized) return;
  myreg<long>();
  myreg<std::complex<double> >();
  myreg<double>();
  initialized = true;
 }

}
}
#endif


