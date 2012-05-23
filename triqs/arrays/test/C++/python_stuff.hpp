
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

#ifndef PYTHON_STUFF_H
#define PYTHON_STUFF_H

#ifndef TRIQS_ARRAYS_WITH_PYTHON_SUPPORT
inline void init_python_stuff(int argc, char **argv) { 
}
#else

#define AS_STRING(X) AS_STRING2(X)
#define AS_STRING2(X) #X

#include <Python.h>
//#include "./src/python/converters.hpp"

//extern "C" { void initmpi();}; 

//void translatorString (std::string const& x) { PyErr_SetString(PyExc_RuntimeError, x.c_str());};

inline void init_python_stuff(int argc, char **argv) { 

 //PyImport_AppendInittab((char*)("mpi"), &initmpi);
 //Py_Initialize(); // init the interpreter
 //triqs::arrays::register_boost_converters();

 // register the exception translation
 //boost::python::register_exception_translator<std::string>(translatorString);

}


#endif
#endif



