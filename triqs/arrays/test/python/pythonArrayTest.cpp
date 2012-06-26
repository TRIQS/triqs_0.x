
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

#include "./test1.cpp"
#include "./src/python/converters.hpp"

//extern "C" { void initmpi();}; 
extern "C" { void init_testarray();}; 

#define AS_STRING(X) AS_STRING2(X)
#define AS_STRING2(X) #X

bool execute_file( const char * filename ) { 
  FILE *fp = fopen(filename,"r");
  if (!fp)  return false;
  PyErr_Clear();
  if (!PyRun_SimpleFile(fp,(char *)filename)) {PyErr_Print(); return false;}
  return true;
 }

int main(int argc, char **argv)
{
 //PyImport_AppendInittab((char*)("mpi"), &initmpi);
 PyImport_AppendInittab((char*)("_testarray"), &init_testarray);
 Py_Initialize(); // init the interpreter
 triqs::arrays::register_boost_converters();
 execute_file ( (std::string(AS_STRING(SOURCE_DIR)) + std::string("/test.py")).c_str() ) ;
 // now start the interpreter main loop and start executing the commands given to the new interpreter...
 //  Py_Main(argc,argv);
}


